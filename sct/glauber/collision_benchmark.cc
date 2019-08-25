#include "sct/glauber/collision.h"
#include "sct/glauber/nucleus.h"
#include "sct/lib/logging.h"

#include "benchmark/benchmark.h"

static void BM_Collision(benchmark::State& state) {
  sct::Nucleus nucleusA;
  sct::parameter_list params;
  params["radius"] = 6.0;
  params["skin_depth"] = 0.5;
  nucleusA.setParameters(state.range(0), params, sct::NucleonPDF::PDF::WoodsSaxon1D);
  nucleusA.generate();
  sct::Nucleus nucleusB;
  nucleusB.setParameters(state.range(0), params, sct::NucleonPDF::PDF::WoodsSaxon1D);
  nucleusB.generate();

  sct::Collision collision;
  collision.setNNCrossSection(4.1);

  for (auto _ : state) {
    collision.collide(nucleusA, nucleusB);
  }
}

static void BM_XSec(benchmark::State& state) {
  sct::parameter_list params;
  params["radius"] = 6.0;
  params["skin_depth"] = 0.5;
  sct::Nucleus nucleusA;
  nucleusA.setParameters(200, params, sct::NucleonPDF::PDF::WoodsSaxon1D);
  nucleusA.generate();
  sct::Nucleus nucleusB;
  nucleusB.setParameters(200, params, sct::NucleonPDF::PDF::WoodsSaxon1D);
  nucleusB.generate();

  sct::Collision collision;
  collision.setNNCrossSection(state.range(0) / 100.0);

  for (auto _ : state) {
    collision.collide(nucleusA, nucleusB);
  }
}

BENCHMARK(BM_Collision)->Range(10, 1000);
BENCHMARK(BM_XSec)->Range(10, 1000);

BENCHMARK_MAIN();
