#include "gtest/gtest.h"
#include "sct/systematics/histogram_collection.h"

#include <string>

#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

using std::string;

// test creation of histograms
TEST(HistogramCollection, creation) {
  sct::HistogramCollection<TH1D> collection;

  collection.add("h", "htest",100, 0, 10);
  collection.add("h2", "hello",20, 1, 100);

  EXPECT_NE(collection.get("h"), nullptr);
  EXPECT_NE(collection.get("h2"), nullptr);
  EXPECT_EQ(collection.get("h3"), nullptr);

  sct::HistogramCollection<TH2D> collection_2d;

  collection_2d.add("h2d", "2d", 100, 0, 1, 100, 1, 10);

  EXPECT_NE(collection_2d.get("h2d"), nullptr);
  EXPECT_EQ(collection_2d.get("h2d2"), nullptr);

  sct::HistogramCollection<TProfile> collection_prof;

  collection_prof.add("hprof", "hproftest",50, 0, 10);

  EXPECT_NE(collection_prof.get("hprof"), nullptr);
  EXPECT_EQ(collection_prof.get("hprof2"), nullptr);
}

TEST(HistogramCollection, initialization) {
  sct::HistogramCollection<TH1D> c;

  c.add("h", "htest",100, 0, 10);
  c.add("h2", "hello",20, 1, 100);

  EXPECT_EQ(string(c.get("h")->GetName()), "h");
  EXPECT_EQ(string(c.get("h2")->GetName()), "h2");
  EXPECT_EQ(string(c.get("h")->GetTitle()), "htest");
  EXPECT_EQ(string(c.get("h2")->GetTitle()), "hello");

  EXPECT_EQ(c.get("h")->GetNbinsX(), 100);
  EXPECT_EQ(c.get("h2")->GetNbinsX(), 20);

  sct::HistogramCollection<TH2D> c_2d;

  c_2d.add("h2d", "2d", 100, 0, 1, 101, 1, 10);  

  EXPECT_EQ(string(c_2d.get("h2d")->GetName()), "h2d");
  EXPECT_EQ(string(c_2d.get("h2d")->GetTitle()), "2d");
  EXPECT_EQ(c_2d.get("h2d")->GetNbinsX(), 100);
  EXPECT_EQ(c_2d.get("h2d")->GetNbinsY(), 101);
  EXPECT_NEAR(c_2d.get("h2d")->GetXaxis()->GetXmax(), 1.0, 1e-10);
  EXPECT_NEAR(c_2d.get("h2d")->GetXaxis()->GetXmin(), 0.0, 1e-10);
  EXPECT_NEAR(c_2d.get("h2d")->GetYaxis()->GetXmax(), 10.0, 1e-10);
  EXPECT_NEAR(c_2d.get("h2d")->GetYaxis()->GetXmin(), 1.0, 1e-10);
}

TEST(HistogramCollection, fill) {
  sct::HistogramCollection<TH1D> c;

  int nbins = 10;
  c.add("h", "htest",nbins, 0, nbins);
  c.fill("h", 0.5);
  c.fill("h", 5.5, 2.0);

  for (int i = 1; i < nbins; ++i) {
    switch (i) {
      case 1 :
        EXPECT_NEAR(c.get("h")->GetBinContent(1), 1.0, 1e-10);
      case 6 :
        EXPECT_NEAR(c.get("h")->GetBinContent(6), 2.0, 1e-10);
    }
  }
}


