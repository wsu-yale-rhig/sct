#include "sct/glauber/mc_glauber.h"

#include "benchmark/benchmark.h"

static void BM_ImpactParameter(benchmark::State &state) {
  sct::MCGlauber generator;
  generator.setImpactParameterRange(0.0, state.range(0));
  for (auto _ : state) {
    generator.run(100);
  }
}

static void BM_NucleusSize(benchmark::State &state) {
  sct::parameter_list params;
  params["radius"] = 6.0;
  params["skin_depth"] = 0.5;
  sct::MCGlauber generator(state.range(0), sct::NucleonPDF::PDF::WoodsSaxon1D,
                           params, 4.1, 200);
  for (auto _ : state) {
    generator.run(100);
  }
}

BENCHMARK(BM_ImpactParameter)->Range(1, 20);
BENCHMARK(BM_NucleusSize)->Range(1, 500);

BENCHMARK_MAIN();
