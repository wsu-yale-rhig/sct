#include "sct/lib/flags.h"

namespace sct {
  
  const char* ProgramUsage() {
    return gflags::ProgramUsage();
  }
  
  void SetUsageMessage(const string& str) {
    gflags::SetUsageMessage(str);
  }
  
  bool ParseCommandLineFlags(int* argc, char** argv) {
    return gflags::ParseCommandLineFlags(argc, &argv, true);
  }
  
} //namespace sct
