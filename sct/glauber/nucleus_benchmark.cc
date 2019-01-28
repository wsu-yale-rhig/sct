#include "sct/glauber/nucleus.h"

#include "benchmark/benchmark.h"

static void BM_NucleusCreation(benchmark::State& state) {
  for (auto _ : state) {
    sct::Nucleus nucleus;
  }
}

static void BM_NucleusGeneration(benchmark::State& state) {
  sct::Nucleus nucleus;
  nucleus.setParameters(state.range(0), 6.0, 0.5, 0.0, 0.0);
  for (auto _ : state) {
    nucleus.generate();
  }
}

static void BM_NucleonRepulsion(benchmark::State& state) {
  sct::Nucleus nucleus;
  nucleus.setParameters(500, 6.0, 0.5, 0.0, 0.0);
  for (auto _ : state) {
    nucleus.setRepulsionDistance((double)state.range(0) / 1000.0);
    nucleus.generate();
  }
}

static void BM_NucleusGenerationRepulsion(benchmark::State& state) {
  sct::Nucleus nucleus;
  nucleus.setParameters(state.range(0), 6.0, 0.5, 0.0, 0.0);
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
