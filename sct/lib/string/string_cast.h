#ifndef SCT_LIB_STRING_STRING_CAST_H
#define SCT_LIB_STRING_STRING_CAST_H

#include <sstream>
#include <string>

namespace sct {

template <typename T> bool CanCast(std::string s) {
  std::istringstream iss(s);
  T dummy;
  iss >> std::skipws >> dummy;
  return iss && iss.eof();
}

template <typename T> T CastTo(std::string s) {
  std::istringstream iss(s);
  T dummy;
  iss >> std::skipws >> dummy;
  return dummy;
}

template <typename T> std::vector<T> ParseStrToVec(std::string str) {
  std::vector<T> ret;
  std::string token;
  while (str.find(" ") != std::string::npos) {
    size_t pos = str.find(" ");
    token = str.substr(0, pos);
    if (CanCast<T>(token)) {
      ret.push_back(CastTo<T>(token));
    }
    str.erase(0, pos + 1);
  }
  if (CanCast<T>(str))
    ret.push_back(CastTo<T>(str));

  return ret;
}
} // namespace sct

#endif // SCT_LIB_STRING_STRING_CAST_H