#include "reader.h"
#include <glog/logging.h>

Reader::Reader(std::istream& input) : input_(input) {
}

bool Reader::ReadBit() {
    if (buffer_index_ == 8) {
        buffer_ = input_.get();
        if (buffer_ == 0xff) {
            uint8_t x = input_.get();
            if (x != 0) {
                LOG(INFO) << int(x);
                throw std::invalid_argument("Wanted to see 0x00 after 0xff");
            }
        }
        buffer_index_ = 0;
    }
    return (buffer_ >> (7 - buffer_index_++)) & 1;
}
