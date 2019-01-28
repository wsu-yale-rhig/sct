#ifndef SCT_LIB_MEMORY_H
#define SCT_LIB_MEMORY_H

// memory includes for libsct's internal use

#include <memory>

// In case we're compiling with c++11, we will emulate c++14's std::make_unique
#if __cplusplus < 201402L && (!defined __cpp_lib_make_unique)
#include "sct/lib/memory/make_unique.h"
#endif

namespace sct {
// commonly used memory functionality
using std::make_shared;
using std::make_unique;
using std::shared_ptr;
using std::unique_ptr;
}  // namespace sct

#endif  // SCT_LIB_MEMORY_H
