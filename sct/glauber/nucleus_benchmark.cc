#include "sct/glauber/nucleus.h"

#include "benchmark/benchmark.h"

static void BM_NucleusCreation(benchmark::State &state) {
  for (auto _ : state) {
    sct::Nucleus nucleus;
  }
}

static void BM_NucleusGeneration(benchmark::State &state) {
  sct::Nucleus nucleus;
  sct::parameter_list params;
  params["radius"] = 6.0;
  params["skin_depth"] = 0.5;
  nucleus.setParameters(state.range(0), params,
                        sct::NucleonPDF::PDF::WoodsSaxon1D);
  for (auto _ : state) {
    nucleus.generate();
  }
}

static void BM_NucleonRepulsion(benchmark::State &state) {
  sct::Nucleus nucleus;
  sct::parameter_list params;
  params["radius"] = 6.0;
  params["skin_depth"] = 0.5;
  nucleus.setParameters(500, params,
                        sct::NucleonPDF::PDF::WoodsSaxon1D);
  for (auto _ : state) {
    nucleus.setRepulsionDistance((double)state.range(0) / 1000.0);
    nucleus.generate();
  }
}

static void BM_NucleusGenerationRepulsion(benchmark::State &state) {
  sct::Nucleus nucleus;
  sct::parameter_list params;
  params["radius"] = 6.0;
  params["skin_depth"] = 0.5;
  nucleus.setParameters(state.range(0), params,
                        sct::NucleonPDF::PDF::WoodsSaxon1D);
  nucleus.setRepulsionDistance((double)state.range(1) / 1000.0);
  for (auto _ : state) {
    nucleus.generate();
  }
}

BENCHMARK(BM_NucleusCreation);
BENCHMARK(BM_NucleusGeneration)->Range(10, 1000);
BENCHMARK(BM_NucleonRepulsion)->Range(10, 1000);
BENCHMARK(BM_NucleusGenerationRepulsion)->Ranges({{10, 1000}, {1, 500}});

BENCHMARK_MAIN();
