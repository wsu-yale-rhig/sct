#ifndef SCT_LIB_LOGGING_H
#define SCT_LIB_LOGGING_H

/* sct uses GLog for logging
 */
#include "glog/logging.h"
#include "glog/stl_logging.h"

#include "sct/lib/flags.h"

// used to control levels of output using glog
// levels: INFO = 0, WARNING = 1, ERROR = 2,
// FATAL = 3. Fatal will cause exit if NDEBUG is turned on
SCT_DECLARE_int(SCT_log_level);
SCT_DECLARE_int(minloglevel);
// controls verbose logging levels
SCT_DECLARE_int(v);
// logs all output to stderr
SCT_DECLARE_bool(logtostderr);

namespace sct {

bool InitLogging(int* argc, char** argv);
bool InitLogging(char* argv);

}  // namespace sct

#endif  // SCT_LIB_LOGGING_H
