#ifndef SCT_LIB_STRING_STRING_CAST_H
#define SCT_LIB_STRING_STRING_CAST_H

#include <sstream>
#include <string>

namespace sct {

template <typename T>
bool CanCast(std::string s) {
  std::istringstream iss(s);
  T dummy;
  iss >> std::skipws >> dummy;
  return iss && iss.eof();
}

template <typename T>
T CastTo(std::string s) {
  std::istringstream iss(s);
  T dummy;
  iss >> std::skipws >> dummy;
  return dummy;
}

}  // namespace sct

#endif  // SCT_LIB_STRING_STRING_CAST_H
