#pragma once

#include <stdint.h>

namespace markers {
const uint16_t kStart = 0xffd8;
const uint16_t kEnd = 0xffd9;
const uint16_t kComment = 0xfffe;
const uint16_t kAppl = 0xffe0;
const uint16_t kAppr = 0xffef;
const uint16_t kDqt = 0xffdb;
const uint16_t kSof0 = 0xffc0;
const uint16_t kDht = 0xffc4;
const uint16_t kSos = 0xffda;
}  // namespace markers
