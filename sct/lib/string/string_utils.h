#ifndef SCT_LIB_STRING_STRING_UTIL_H
#define SCT_LIB_STRING_STRING_UTIL_H

// common string functionality

#include <sstream>
#include <set>
#include <vector>

#include "sct/lib/string/string_cast.h"
#include "sct/lib/string/string.h"

namespace sct {
  
  // allows the user to split a string on a specified delimiter
  // default delimiter is ' '
  
  template <class Container>
  void SplitString(const string& str, Container& cont, char delim = ' ') {
    std::stringstream ss(str);
    string token;
    while (std::getline(ss, token, delim)) {
      cont.push_back(token);
    }
  }
  
  template<>
  void SplitString(const string& str, std::set<string>& cont, char delim);
  
  // template function to make a string out of an arbitrary number
  // of varied inputs, as long as the input has an sstream << overload
  // for instance MakeString("string", 42, 36.5, ClassWith_<<_Overload);
  
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
  
  // remove all instances of pred in a string
  
  template<typename T, typename P>
  T RemoveIf(T beg, T end, P pred)
  {
    T dest = beg;
    for (T itr = beg;itr != end; ++itr)
      if (!pred(*itr))
        *(dest++) = *itr;
    return dest;
  }
  
  // checks if full_string ends with ending
  bool HasEnding(const string& full_string, const string& ending);
  
  // if a string begins with seq, removes seq in-place and returns true
  // otherwise returns false and does not modify arg
  bool Consume(string& arg, const string& seq);
  
  // if a string begins with c, removes c in-place and returns true
  // otherwise returns false and does not modify arg
  bool Consume(string& arg, char c);
  
  // returns a string containing every character up till the next occurence of sec
  string SplitOnNextOccurence(string& arg, const string& seq);
  
  // parses a comma separated string and casts to type T,
  // returning a set
  template<typename T, class Container=std::set<T>>
  Container ParseArgString(string str) {
    Container ret;
    
    // remove all spaces
    str.erase(::sct::RemoveIf(str.begin(), str.end(), ::isspace), str.end());
    
    string token;
    while ( str.find(",") != string::npos ) {
      size_t pos = str.find(",");
      token = str.substr(0, pos);
      if (CanCast<T>(token)) {
        ret.insert(CastTo<T>(token));
        str.erase(0, pos + 1);
      }
    }
    if (CanCast<T>(str))
      ret.insert(CastTo<T>(str));
    
    return ret;
  }
  
  template<typename T>
  std::vector<T> ParseArgStringToVec(string str) {
    std::vector<T> ret;
    
    // remove all spaces
    str.erase(::sct::RemoveIf(str.begin(), str.end(), ::isspace), str.end());
    
    string token;
    while ( str.find(",") != string::npos ) {
      size_t pos = str.find(",");
      token = str.substr(0, pos);
      if (CanCast<T>(token)) {
        ret.push_back(CastTo<T>(token));
        str.erase(0, pos + 1);
      }
    }
    if (CanCast<T>(str))
      ret.push_back(CastTo<T>(str));
    
    return ret;
  }

  template <class Container>
  string Concat(const Container& c, const string& delim) {
    std::stringstream ss;
    int count = c.size() - 1;
    for (auto i = c.begin(); i != c.end(); ++i, --count) {
      ss << *i << (count ? delim : "");
    }
    return string(ss.str());
  }
  
}// namespace sct

#endif // SCT_LIB_STRING_STRING_UTIL_H

