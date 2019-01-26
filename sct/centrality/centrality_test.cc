#include "gtest/gtest.h"
#include "sct/centrality/centrality.h"

#include <random>

// first test that the centrality integration works
TEST(Centrality, centralityIntegration0_100) {
  
  std::unique_ptr<TH1D> h = sct::make_unique<TH1D>("h", "", 150, 0, 150);
  
  std::random_device random;
  std::mt19937 generator(random());
  std::uniform_int_distribution<> dis(0, 100);
  for (int i = 1; i <= 1e7; ++i) {
    h->Fill(dis(generator));
  }
  
  sct::Centrality centrality;
  centrality.setSimuRefmult(h.get());
  std::vector<unsigned> boundaries_forward = centrality.centralityBins();
  std::vector<unsigned> boundaries_backward = centrality.centralityBins();
  
  EXPECT_EQ(boundaries_forward, boundaries_backward);
  EXPECT_EQ(boundaries_forward, std::vector<unsigned>({20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95}));
  
}

// same test, but for a random distribution
TEST(Centrality, centralityIntegrationRandomDist) {
  
  std::unique_ptr<TH1D> h = sct::make_unique<TH1D>("h", "", 150, 0, 150);
  std::random_device random;
  std::mt19937 generator(random());
  std::uniform_int_distribution<> dis(0, 150);
  for (int i = 1; i <= 1000; ++i) {
    h->Fill(dis(generator));
  }
  
  sct::Centrality centrality;
  centrality.setSimuRefmult(h.get());
  std::vector<unsigned> boundaries_forward = centrality.centralityBins();
  std::vector<unsigned> boundaries_backward = centrality.centralityBins();
  
  EXPECT_EQ(boundaries_forward, boundaries_backward);
  
}
