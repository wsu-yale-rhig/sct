#ifndef SCT_LIB_STRING_PATH_UTILS_H
#define SCT_LIB_STRING_PATH_UTILS_H

// functionality to manipulate strings representing paths
// and file names

#include "sct/lib/string/string.h"

namespace sct {

// strips the path and returns the file name
string GetFileName(const string& path);

// strips the file name and returns the path
string GetPath(const string& path);

}  // namespace sct

#endif  // SCT_LIB_STRING_PATH_UTILS_H
