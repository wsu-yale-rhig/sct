#ifndef SCT_LIB_FLAGS_H
#define SCT_LIB_FLAGS_H

/* implementation of command line argument parser
 * using gflags
 */

#include "sct/lib/string/string.h"

namespace sct {

// returns the usage message
const char* ProgramUsage();

// sets the usage message that is printed when the executable
// is called with flag --help. Do not add lines for individual
// flags, this is handled by gflags
void SetUsageMessage(const string& str);

// Parses command line options. Follows GFlags convention,
// where argv and argc are modified to remove any flags that
// have been handled, leaving any unhandled flags
bool ParseCommandLineFlags(int* argc, char** argv);

}  // namespace sct

#include "gflags/gflags.h"

#define SCT_GFLAGS_DEFINE(type, name, default_value, help_str) \
  DEFINE_##type(name, default_value, help_str);                \
  namespace sct {                                              \
  using ::FLAGS_##name;                                        \
  }  // namespace sct

#define SCT_DEFINE_int(name, default_value, help_str) \
  SCT_GFLAGS_DEFINE(int32, name, default_value, help_str)
#define SCT_DEFINE_int64(name, default_value, help_str) \
  SCT_GFLAGS_DEFINE(int64, name, default_value, help_str)
#define SCT_DEFINE_double(name, default_value, help_str) \
  SCT_GFLAGS_DEFINE(double, name, default_value, help_str)
#define SCT_DEFINE_bool(name, default_value, help_str) \
  SCT_GFLAGS_DEFINE(bool, name, default_value, help_str)
#define SCT_DEFINE_string(name, default_value, help_str) \
  SCT_GFLAGS_DEFINE(string, name, default_value, help_str)

#define SCT_GFLAGS_DECLARE(type, name) \
  DECLARE_##type(name);                \
  namespace sct {                      \
  using ::FLAGS_##name;                \
  }  // namespace sct

#define SCT_DECLARE_int(name) SCT_GFLAGS_DECLARE(int32, name)
#define SCT_DECLARE_int64(name) SCT_GFLAGS_DECLARE(int64, name)
#define SCT_DECLARE_double(name) SCT_GFLAGS_DECLARE(double, name)
#define SCT_DECLARE_bool(name) SCT_GFLAGS_DECLARE(bool, name)
#define SCT_DECLARE_string(name) SCT_GFLAGS_DECLARE(string, name)

#endif  // SCT_LIB_FLAGS_H
