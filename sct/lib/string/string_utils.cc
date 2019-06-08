#include "sct/lib/string/string_utils.h"

namespace sct {

template <>
void SplitString(const string &str, std::set<string> &cont, char delim) {
  std::stringstream ss(str);
  string token;
  while (std::getline(ss, token, delim)) {
    cont.insert(token);
  }
}

template <> string MakeString(const string &s) { return s; }

string MakeString(const char *chr) { return std::string(chr); }

bool HasEnding(const string &full_string, const string &ending) {
  if (full_string.length() >= ending.length()) {
    return (0 == full_string.compare(full_string.length() - ending.length(),
                                     ending.length(), ending));
  } else {
    return false;
  }
}

bool Consume(string &arg, const string &seq) {
  if (seq.size() > arg.size())
    return false;
  if (arg.compare(0, seq.length(), seq) == 0) {
    arg = string(arg.begin() + seq.length(), arg.end());
    return true;
  }
  return false;
}

bool Consume(string &arg, char c) {
  if (arg.size() == 0)
    return false;
  if (arg[0] == c) {
    arg = string(arg.begin() + 1, arg.end());
    return true;
  }
  return false;
}

string SplitOnNextOccurence(string &arg, const string &seq) {
  size_t loc = arg.find(seq);
  if (loc == string::npos)
    return string();
  string ret = string(arg.begin(), arg.begin() + loc);
  arg = string(arg.begin() + loc + seq.size(), arg.end());
  return ret;
}

} // namespace sct
