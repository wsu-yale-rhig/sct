// sct/systematics/histogram_collection.hh

#ifndef SCT_SYSTEMATICS_HISTOGRAM_COLLECTION_HH
#define SCT_SYSTEMATICS_HISTOGRAM_COLLECTION_HH

#include <string>
#include <unordered_map>
#include <memory>
#include <iostream>

#include "sct/core/base.hh"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TProfile2D.h"

namespace sct {
    
  template<class H, class Key=std::string, class hash=std::hash<Key>>
  class HistogramCollection {
  public:
    HistogramCollection() : histograms_() { };
    ~HistogramCollection() { };
    
    H* get(Key key) {
      if (keyExists(key))
        return histograms_[key].get();
      return nullptr;
    }
    
    template <typename... Args>
    void add(Key key, Args... args) {
      histograms_[key] = make_unique<H>(key.c_str(), args...);
      histograms_[key]->SetDirectory(0);
    }
    
    template <typename ...Args>
    bool fill(Key key, Args... args) {
      if (!keyExists(key))
        return false;
      histograms_[key]->Fill(args...);
      return true;
    }
    
    void write() {
      for (auto& h : histograms_)
        h.second->Write();
    }
    
    void clear() {
      histograms_.clear();
    }
    
  private:
    
    unordered_map<Key, unique_ptr<H>, hash> histograms_ = {};
    
    bool keyExists(string key) {
      for (auto& h : histograms_) {
        if (h.first == key)
        return true;
      }
      return false;
    }
    
  };
  
} // namespace sct

#endif // SCT_SYSTEMATICS_HISTOGRAM_COLLECTION_HH
