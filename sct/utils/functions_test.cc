#include "gtest/gtest.h"

#include "sct/lib/logging.h"
#include "sct/lib/math.h"
#include "sct/lib/memory.h"
#include "sct/utils/functions.h"

#include <vector>

#include "TF1.h"
#include "TF3.h"

TEST(function, step_function_1d) {

  sct::unique_ptr<TF1> f1 =
      sct::make_unique<TF1>("testf", sct::StepFunction1D, 1, 0, 5);
  f1->SetParameter(0, 1.0);

  EXPECT_NEAR(f1->Eval(0.5), 1.0, 1e-5);
  EXPECT_NEAR(f1->Eval(0.9999999), 1.0, 1e-5);
  EXPECT_NEAR(f1->Eval(1.5), 0.0, 1e-5);
  EXPECT_NEAR(f1->Eval(1.0000001), 0.0, 1e-5);
}

TEST(function, step_function_3d) {

  sct::unique_ptr<TF3> f1 =
      sct::make_unique<TF3>("testf", sct::StepFunction, 1, -5, 5, -5, 5, -5, 5);

  double max_r = 1.5;
  double area = sct::pi * pow(max_r, 2.0);
  double volume = 4.0 / 3.0 * sct::pi * pow(max_r, 3.0);
  f1->SetParameter(0, area);

  std::vector<std::vector<double>> in{{1.4999, 0.0, 0.0},
                                      {0.0, 1.49999, 0.0},
                                      {0.0, 0.0, 1.4999},
                                      {1.4999, 0.000001, 0.0},
                                      {1.0, 1.0, 0.0},
                                      {0.74, 0.74, 0.74}};
  std::vector<std::vector<double>> out{{1.50001, 0.0, 0.0},
                                       {0.0, 1.50001, 0.0},
                                       {0.0, 0.0, 1.50001},
                                       {1.4999, 0.5, 0.0}};

  for (auto& e : in) {
    EXPECT_NEAR(f1->Eval(e[0], e[1], e[2]), 1.0/volume, 1e-5);
  }
  for (auto& e : out) {
    EXPECT_NEAR(f1->Eval(e[0], e[1], e[2]), 0.0, 1e-5);
  }
}