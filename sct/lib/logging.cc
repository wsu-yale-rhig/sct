#include "sct/lib/logging.h"

// initial setting for logging level
SCT_DEFINE_int(sct_log_level, google::ERROR,
               "Minimum severity level to output.");

namespace sct {

bool InitLogging(int* argc, char** argv) {
  if (*argc == 0) return true;
  ::google::InitGoogleLogging(argv[0]);
  ::google::InstallFailureSignalHandler();

  // choose the lower of sct log level & minimum log level
  FLAGS_minloglevel = std::min(FLAGS_sct_log_level, FLAGS_minloglevel);

  return true;
}

bool InitLogging(char* argv) {
  if (argv == nullptr) return true;
  ::google::InitGoogleLogging(argv);
  ::google::InstallFailureSignalHandler();

  // choose the lower of sct log level & minimum log level
  FLAGS_minloglevel = std::min(FLAGS_sct_log_level, FLAGS_minloglevel);

  return true;
}

}  // namespace sct
