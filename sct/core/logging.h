#ifndef SCT_CORE_LOGGING_H
#define SCT_CORE_LOGGING_H

/* sct uses GLog for logging
 */
#include "glog/stl_logging.h"
#include "glog/logging.h"

#include "sct/core/flags.h"


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
  
  // string utilities
  template <class Container>
  string Concat(const Container& c, const string& delim) {
    std::stringstream ss;
    int count = c.size() - 1;
    for (auto i = c.begin(); i != c.end(); ++i, --count) {
      ss << *i << (count ? delim : "");
    }
    return string(ss.str());
  }
  
  inline void MakeStringInternal(std::stringstream& ss) {
    return;
  }
  
  template <typename T>
  inline void MakeStringInternal(std::stringstream& ss, T& t) {
    ss << t;
  }
  
  template <typename T, typename... Args>
  inline void MakeStringInternal(std::stringstream& ss, T& t, const Args&... args) {
    MakeStringInternal(ss, t);
    MakeStringInternal(ss, args...);
  }
  
  template <typename... Args>
  string MakeString(const Args&... args) {
    std::stringstream ss;
    MakeStringInternal(ss, args...);
    return string(ss.str());
  }
  
  template<>
  string MakeString(const string& s);
  
  string MakeString(const char* chr);
  
  template <class Container>
  void SplitString(const std::string& str, Container& cont, char delim = ' ') {
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delim)) {
      cont.push_back(token);
    }
  }

  
} // namespace sct

#endif // SCT_CORE_LOGGING_H
