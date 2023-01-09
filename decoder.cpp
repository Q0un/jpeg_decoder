#include <decoder.h>
#include <glog/logging.h>

#include <array>
#include <vector>
#include <unordered_set>
#include <cmath>
#include <limits>

#include "huffman.h"
#include "reader.h"
#include "markers.h"
#include "fft.h"
#include "zigzag.h"

namespace {

void ReadComment(Reader& reader, Image& image, bool need = false) {
    uint16_t len = reader.ReadWord<2>();
    std::string comment;
    for (uint16_t i = 0; i < len - 2; ++i) {
        comment += reader.ReadWord<1>();
    }
    if (need) {
        image.SetComment(comment);
    }
}

using QuantMatrix = std::array<std::array<uint16_t, 8>, 8>;

void ReadQuantTables(Reader& reader, std::array<QuantMatrix, 16>& tables) {
    uint16_t len = reader.ReadWord<2>();
    len -= 2;
    while (len) {
        uint8_t info = reader.ReadWord<1>();
        uint8_t id = info & ((1 << 4) - 1);
        uint8_t value_len = (info ^ id) >> 4;
        if (value_len == 0) {
            len -= 65;
        } else if (value_len == 1) {
            len -= 129;
        } else {
            throw std::invalid_argument("Wrong file format");
        }
        for (const auto& [i, j] : kZigzag) {
            if (value_len == 0) {
                tables[id][i][j] = reader.ReadWord<1>();
            } else {
                tables[id][i][j] = reader.ReadWord<2>();
            }
        }
    }
}

struct ChannelInfo {
    uint8_t id = 255;
    uint8_t thinning;
    uint8_t quant_table_id;
    uint8_t ac_id;
    uint8_t dc_id;
};

struct MetaInfo {
    uint8_t precision;
    uint16_t height;
    uint16_t width;
    uint8_t channels_count;
    std::vector<ChannelInfo> channels;
};

MetaInfo ReadMetaInfo(Reader& reader) {
    MetaInfo info;
    uint16_t len = reader.ReadWord<2>();
    info.precision = reader.ReadWord<1>();
    info.height = reader.ReadWord<2>();
    info.width = reader.ReadWord<2>();
    uint8_t channels_count = reader.ReadWord<1>();
    info.channels_count = channels_count;
    info.channels.resize(3);
    for (uint8_t i = 0; i < channels_count; ++i) {
        ChannelInfo channel;
        channel.id = reader.ReadWord<1>() - 1;
        channel.thinning = reader.ReadWord<1>();
        channel.quant_table_id = reader.ReadWord<1>();
        if (channel.id != 0 && channel.id != 1 && channel.id != 2) {
            throw std::invalid_argument("Wrong channel id");
        }
        if (channel.quant_table_id >= 16) {
            throw std::invalid_argument("Wrong quant table id");
        }
        info.channels[channel.id] = std::move(channel);
    }
    return info;
}

struct HuffmanTable {
    std::vector<uint8_t> code_lengths;
    std::vector<uint8_t> values;
};

void ReadHuffmanTables(Reader& reader, std::array<std::array<HuffmanTree, 2>, 16>& huffmans) {
    uint16_t len = reader.ReadWord<2>();
    len -= 2;
    while (len) {
        uint8_t id = reader.ReadWord<1>();
        --len;
        HuffmanTable table;
        table.code_lengths.resize(16);
        size_t count_values = 0;
        for (size_t i = 0; i < 16; ++i) {
            table.code_lengths[i] = reader.ReadWord<1>();
            count_values += table.code_lengths[i];
            --len;
        }
        for (uint16_t i = 0; i < count_values; ++i) {
            table.values.push_back(reader.ReadWord<1>());
            --len;
        }
        HuffmanTree tree;
        tree.Build(table.code_lengths, table.values);
        uint8_t buf = id & ((1 << 4) - 1);
        uint8_t huf_class = (id ^ buf) >> 4;
        if (huf_class > 1) {
            throw std::invalid_argument("Wrong huffman table class");
        }
        id = buf;
        huffmans[id][huf_class] = std::move(tree);
    }
}

void ReadSOSInfo(Reader& reader, MetaInfo& meta_info) {
    uint16_t len = reader.ReadWord<2>();
    uint8_t channels = reader.ReadWord<1>();
    if (channels != meta_info.channels_count) {
        throw std::invalid_argument("Wrong number of channels");
    }
    for (uint8_t i = 0; i < channels; ++i) {
        uint8_t ch_id = reader.ReadWord<1>();
        uint8_t ids = reader.ReadWord<1>();
        bool fl = false;
        for (const auto& ch : meta_info.channels) {
            if (ch.id + 1 == ch_id) {
                fl = true;
            }
        }
        if (!fl) {
            throw std::invalid_argument("Wrong channel id");
        }
        meta_info.channels[ch_id - 1].ac_id = ids & ((1 << 4) - 1);
        meta_info.channels[ch_id - 1].dc_id = (ids ^ meta_info.channels[ch_id - 1].ac_id) >> 4;
        if (meta_info.channels[ch_id - 1].ac_id >= 16 ||
            meta_info.channels[ch_id - 1].dc_id >= 16) {
            throw std::invalid_argument("Wrong huffman table id");
        }
    }
    if (len - 3 - channels * 2 != 3) {
        throw std::invalid_argument("Wrong progressive stuff");
    }
    if (reader.ReadWord<1>() != 0) {
        throw std::invalid_argument("Wrong progressive stuff");
    }
    if (reader.ReadWord<1>() != 63) {
        throw std::invalid_argument("Wrong progressive stuff");
    }
    if (reader.ReadWord<1>() != 0) {
        throw std::invalid_argument("Wrong progressive stuff");
    }
}

using DataMatrix = std::array<std::array<int16_t, 8>, 8>;

void ReadDataMatrix(Reader& reader, std::pair<HuffmanTree&, HuffmanTree&> huffmans,
                    DataMatrix& matrix) {
    for (size_t i = 0; i < 8; ++i) {
        for (size_t j = 0; j < 8; ++j) {
            matrix[i][j] = 0;
        }
    }
    for (size_t it = 0; it < 64; ++it) {
        int len;
        uint8_t i = kZigzag[it].first;
        uint8_t j = kZigzag[it].second;
        bool bit;
        if (i == 0 && j == 0) {
            do {
                bit = reader.ReadBit();
            } while (!huffmans.first.Move(bit, len));
            if (len == 0) {
                continue;
            }
            int16_t dc;
            dc = reader.ReadBit();
            bool neg = 0;
            if (dc == 0) {
                neg = 1;
                dc = 1;
            }
            for (uint8_t k = 1; k < len; ++k) {
                dc <<= 1;
                dc |= neg ^ reader.ReadBit();
            }
            if (neg) {
                dc *= -1;
            }
            matrix[0][0] = dc;
        } else {
            do {
                bit = reader.ReadBit();
            } while (!huffmans.second.Move(bit, len));
            if (len == 0) {
                return;
            }
            uint8_t buf = len & ((1 << 4) - 1);
            uint8_t zeros = (len ^ buf) >> 4;
            len = buf;
            for (uint8_t k = 0; k < zeros; ++k) {
                matrix[i][j] = 0;
                ++it;
                if (it == 64) {
                    throw std::invalid_argument("Wrong data matrix");
                }
                i = kZigzag[it].first;
                j = kZigzag[it].second;
            }
            if (len == 0) {
                continue;
            }
            int16_t ac;
            ac = reader.ReadBit();
            bool neg = 0;
            if (ac == 0) {
                neg = 1;
                ac = 1;
            }
            for (uint8_t k = 1; k < len; ++k) {
                ac <<= 1;
                ac |= neg ^ reader.ReadBit();
            }
            if (neg) {
                ac *= -1;
            }
            matrix[i][j] = ac;
        }
    }
}

int16_t MakeChannelMatrixGood(std::vector<DataMatrix>& ch_datas, QuantMatrix& quant, int16_t last,
                              DctCalculator& dct, std::vector<double>& bad_matrix,
                              std::vector<double>& good_matrix) {

    for (size_t k = 0; k < ch_datas.size(); ++k) {
        ch_datas[k][0][0] += last;
        last = ch_datas[k][0][0];
        for (uint8_t i = 0; i < 8; ++i) {
            for (uint8_t j = 0; j < 8; ++j) {
                ch_datas[k][i][j] *= quant[i][j];
                bad_matrix[i * 8 + j] = ch_datas[k][i][j];
            }
        }

        dct.Inverse();
        for (uint8_t i = 0; i < 8; ++i) {
            for (uint8_t j = 0; j < 8; ++j) {
                if (good_matrix[i * 8 + j] > std::numeric_limits<int16_t>::max() ||
                    good_matrix[i * 8 + j] < std::numeric_limits<int16_t>::lowest()) {
                    throw std::invalid_argument("Numbers are too big");
                }
                ch_datas[k][i][j] = good_matrix[i * 8 + j];
            }
        }
    }
    return last;
}

RGB YCbCrToRGB(int16_t y, int16_t cb, int16_t cr) {
    RGB rgb;
    rgb.r = std::round(y + 1.402 * (cr - 128));
    rgb.g = std::round(y - 0.34414 * (cb - 128) - 0.71414 * (cr - 128));
    rgb.b = std::round(y + 1.772 * (cb - 128));
    rgb.r = std::min(std::max(0, rgb.r), 255);
    rgb.g = std::min(std::max(0, rgb.g), 255);
    rgb.b = std::min(std::max(0, rgb.b), 255);
    return rgb;
}

void ReadData(Reader& reader, std::array<std::array<HuffmanTree, 2>, 16>& huffmans, MetaInfo& info,
              Image& image, std::array<QuantMatrix, 16>& quants) {
    uint8_t hor = 1, vert = 1;
    for (const auto& ch : info.channels) {
        uint8_t local_vert = ch.thinning & ((1 << 4) - 1);
        vert = std::max(vert, local_vert);
        uint8_t local_hor = (ch.thinning ^ local_vert) >> 4;
        hor = std::max(hor, local_hor);
    }
    std::vector<double> bad_matrix(64);
    std::vector<double> good_matrix(64);
    DctCalculator dct(8, &bad_matrix, &good_matrix);
    const size_t h = 8 * vert;
    const size_t w = 8 * hor;
    std::vector<int16_t> lasts(info.channels.size(), 0);
    std::vector<std::vector<std::vector<int16_t>>> datas;
    datas.resize(info.channels_count);
    std::vector<DataMatrix> ch_datas;
    std::vector<std::vector<int16_t>> local_res;
    for (size_t i = 0; i < (info.height + h - 1) / h; ++i) {
        for (size_t j = 0; j < (info.width + w - 1) / w; ++j) {
            for (const auto& ch : info.channels) {
                if (ch.id == 255) {
                    continue;
                }
                uint8_t local_vert = ch.thinning & ((1 << 4) - 1);
                uint8_t local_hor = (ch.thinning ^ local_vert) >> 4;
                if (local_vert == 0 || local_hor == 0) {
                    throw std::invalid_argument("Wrong thinning param");
                }
                if (vert % local_vert != 0 || hor % local_hor != 0) {
                    throw std::invalid_argument("Wrong thinning param");
                }
                ch_datas.resize(local_vert * local_hor);
                for (uint8_t k = 0; k < local_vert * local_hor; ++k) {
                    ReadDataMatrix(reader, {huffmans[ch.dc_id][0], huffmans[ch.ac_id][1]},
                                   ch_datas[k]);
                }
                lasts[ch.id] = MakeChannelMatrixGood(ch_datas, quants[ch.quant_table_id],
                                                     lasts[ch.id], dct, bad_matrix, good_matrix);
                local_res.resize(8 * local_vert);
                for (size_t x = 0; x < 8 * local_vert; ++x) {
                    local_res[x].resize(8 * local_hor);
                    for (size_t y = 0; y < 8 * local_hor; ++y) {
                        local_res[x][y] = ch_datas[(x / 8) * local_hor + (y / 8)][x % 8][y % 8];
                    }
                }
                datas[ch.id] = local_res;
            }
            for (size_t x = 0; x < h; ++x) {
                for (size_t y = 0; y < w; ++y) {
                    uint8_t local_vert = info.channels[0].thinning & ((1 << 4) - 1);
                    uint8_t local_hor = (info.channels[0].thinning ^ local_vert) >> 4;
                    int16_t yel = info.channels[0].id == 255
                                      ? 0
                                      : datas[0][x / (vert / local_vert)][y / (hor / local_hor)];
                    local_vert = info.channels[1].thinning & ((1 << 4) - 1);
                    local_hor = (info.channels[1].thinning ^ local_vert) >> 4;
                    int16_t cb = info.channels[1].id == 255
                                     ? 0
                                     : datas[1][x / (vert / local_vert)][y / (hor / local_hor)];
                    local_vert = info.channels[2].thinning & ((1 << 4) - 1);
                    local_hor = (info.channels[2].thinning ^ local_vert) >> 4;
                    int16_t cr = info.channels[2].id == 255
                                     ? 0
                                     : datas[2][x / (vert / local_vert)][y / (hor / local_hor)];
                    yel = std::min(std::max(0, yel + 128), 255);
                    cb = std::min(std::max(0, cb + 128), 255);
                    cr = std::min(std::max(0, cr + 128), 255);
                    if (i * h + x < image.Height() && j * w + y < image.Width()) {
                        image.SetPixel(i * h + x, j * w + y, YCbCrToRGB(yel, cb, cr));
                    }
                }
            }
        }
    }
}

}  // anonymous namespace

Image Decode(std::istream& input) {
    Reader reader(input);
    uint16_t marker = reader.ReadWord<2>();
    if (input.eof() || marker != markers::kStart) {
        throw std::invalid_argument("Wrong file format");
    }
    Image image;
    std::array<QuantMatrix, 16> quant_matrixes;
    MetaInfo meta_info;
    std::array<std::array<HuffmanTree, 2>, 16> huffmans;
    std::unordered_set<uint16_t> used_markers;
    used_markers.insert(marker);
    while (true) {
        marker = reader.ReadWord<2>();
        if (input.eof()) {
            throw std::invalid_argument("Wrong file format");
        }
        if (marker >= markers::kAppl && marker <= markers::kAppr) {
            ReadComment(reader, image);
            continue;
        }
        if (used_markers.contains(marker) && marker != markers::kDht && marker != markers::kDqt) {
            throw std::invalid_argument("Marker could only be used one time");
        }
        used_markers.insert(marker);
        switch (marker) {
            case markers::kEnd: {
                LOG(INFO) << "Decoding end";
                return image;
            }
            case markers::kComment: {
                ReadComment(reader, image, true);
                LOG(INFO) << "Comment read";
                break;
            }
            case markers::kDqt: {
                ReadQuantTables(reader, quant_matrixes);
                LOG(INFO) << "Quant matrix read";
                break;
            }
            case markers::kSof0: {
                meta_info = ReadMetaInfo(reader);
                image.SetSize(meta_info.width, meta_info.height);
                LOG(INFO) << "Meta info read";
                break;
            }
            case markers::kDht: {
                ReadHuffmanTables(reader, huffmans);
                LOG(INFO) << "Huffman read";
                break;
            }
            case markers::kSos: {
                ReadSOSInfo(reader, meta_info);
                LOG(INFO) << "SOS read";
                ReadData(reader, huffmans, meta_info, image, quant_matrixes);
                LOG(INFO) << "Data read";
                break;
            }
            default:
                throw std::invalid_argument("Wrong file format");
        }
    }
}
