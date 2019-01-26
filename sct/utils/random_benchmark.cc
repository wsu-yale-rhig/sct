#include "benchmark/benchmark.h"
#include "sct/core/logging.h"
#include "sct/utils/random.h"

#include "TH1.h"
#include "TF1.h"

static void BM_random_sample(benchmark::State& state) {
  
  sct::Random& rand = sct::Random::instance();
  double total = 0.0;
  for (auto _ : state) {
    benchmark::DoNotOptimize(total += rand.linear());
  }
}

static void BM_root_th1_sample(benchmark::State& state) {
  
  TH1D* h = new TH1D("h", "h", 100, 0, 100);
  for (int i = 1; i <= 100; ++i)
    h->SetBinContent(i, i);
  
  double total = 0.0;
  for (auto _ : state) {
    benchmark::DoNotOptimize(total += h->GetRandom());
  }
  delete h;
}

static void BM_root_tf1_sample(benchmark::State& state) {
  
  TF1* f = new TF1("f", "[0]*x", 0, 100);
  f->SetParameter(0, 1);
  
  double total = 0.0;
  for (auto _ : state) {
    benchmark::DoNotOptimize(total += f->GetRandom());
  }
  delete f;
}

BENCHMARK(BM_random_sample)->Range(1e4, 1e6);
BENCHMARK(BM_root_th1_sample)->Range(1e4, 1e6);
BENCHMARK(BM_root_tf1_sample)->Range(1e4, 1e6);

BENCHMARK_MAIN();
