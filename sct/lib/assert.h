#ifndef SCT_LIB_ASSERT_H
#define SCT_LIB_ASSERT_H

#include "sct/lib/string/string.h"

#include <exception>
#include <vector>

// verbose exception handling for sct

namespace sct {
// assertion can be anything that will be evaluated to a boolean
// failure will throw an AssertionFailure, which provides some verbose
// error logging (where it was thrown from in the source code, for instance)
#define SCT_ASSERT(assertion, ...)                                 \
  if (!(assertion)) {                                              \
    throw ::sct::AssertionFailure(__FILE__, __LINE__, #assertion,  \
                                  ::sct::MakeString(__VA_ARGS__)); \
  }

// throws an AssertionFailure, as described above
#define SCT_THROW(...)                                             \
  {                                                                \
    throw ::sct::AssertionFailure(__FILE__, __LINE__, "",          \
                                  ::sct::MakeString(__VA_ARGS__)); \
  }

// can collect verbose messages recursively through a call stack
// of try/catch blocks, using Append
class AssertionFailure : public std::exception {
 public:
  AssertionFailure(const char* file, int line, const char* failure,
                   const string& msg);

  const std::vector<string>& MessageStack() const { return msg_stack_; }

  void Append(const string& str);

  string Msg();

  const char* what() const noexcept override;

 private:
  std::vector<string> msg_stack_;
  string msg_;
};

}  // namespace sct

#endif  // SCT_LIB_ASSERT_H