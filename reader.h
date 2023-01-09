#pragma once

#include <stdint.h>
#include <type_traits>
#include <istream>

template <uint8_t N, std::enable_if_t<N == 1 || N == 2>* = nullptr>
using ByteWord = std::conditional_t<N == 1, uint8_t, uint16_t>;

class Reader {
public:
    Reader(std::istream& input);

    template <int N>
    ByteWord<N> ReadWord() {
        ByteWord<N> word;
        input_.read(reinterpret_cast<char*>(&word), sizeof(ByteWord<N>));
        if constexpr (N == 1) {
            return word;
        } else {
            uint16_t buf = word & ((1 << 8) - 1);
            word >>= 8;
            word |= buf << 8 * (N - 1);
            return word;
        }
    }

    bool ReadBit();

private:
    std::istream& input_;
    uint8_t buffer_;
    uint8_t buffer_index_ = 8;
};
