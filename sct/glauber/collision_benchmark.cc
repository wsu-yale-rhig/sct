#include "sct/glauber/collision.h"
#include "sct/glauber/nucleus.h"
#include "sct/lib/logging.h"

#include "benchmark/benchmark.h"

static void BM_Collision(benchmark::State& state) {
  sct::Nucleus nucleusA;
  nucleusA.setParameters(state.range(0), 6.0, 1.0, 0.0, 0.0);
  nucleusA.generate();
  sct::Nucleus nucleusB;
  nucleusB.setParameters(state.range(0), 6.0, 1.0, 0.0, 0.0);
  nucleusB.generate();

  sct::Collision collision;
  collision.setNNCrossSection(4.1);

  for (auto _ : state) {
    collision.collide(nucleusA, nucleusB);
  }
}

static void BM_XSec(benchmark::State& state) {
  sct::Nucleus nucleusA;
  nucleusA.setParameters(200, 6.0, 1.0, 0.0, 0.0);
  nucleusA.generate();
  sct::Nucleus nucleusB;
  nucleusB.setParameters(200, 6.0, 1.0, 0.0, 0.0);
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
