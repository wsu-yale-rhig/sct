#include "sct/lib/assert.h"

#include <algorithm>
#include <numeric>
#include <cstring>

#include "sct/lib/string/string_utils.h"
#include "sct/lib/string/path_utils.h"

namespace sct {
  
  AssertionFailure::AssertionFailure(const char* file,
                                     int line,
                                     const char* failure,
                                     const string& msg) :
    msg_stack_({MakeString("[assertion failure: ", GetFileName(string(file)), "::", line, "] "),
               string(failure), " ", msg}), msg_(Msg()) {}
  
  string AssertionFailure::Msg() {
    return std::accumulate(msg_stack_.begin(), msg_stack_.end(), string(" "));
  }
  
  void AssertionFailure::Append(const string& str) {
    msg_stack_.push_back(str);
    msg_ = Msg();
  }
  
  const char* AssertionFailure::what() const noexcept {
    return msg_.c_str();
  }
  
} // namespace sct
