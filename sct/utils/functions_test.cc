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

TEST(function, step_function_1d_params) {
  EXPECT_EQ(sct::StepFunction1D_npar, 1);
  std::pair<unsigned, double> true_val = {0, 1.0};
  std::pair<unsigned, double> lookup_val = sct::StepFunction1D_params.at("d");
  EXPECT_EQ(true_val, lookup_val);
}

TEST(function, step_function_3d) {

  sct::unique_ptr<TF3> f1 =
      sct::make_unique<TF3>("testf", sct::StepFunction, 1, -5, 5, -5, 5, -5, 5);

  double max_r = 1.5;
  double area = sct::pi * pow(max_r, 2.0);
  double volume = 4.0 / 3.0 * sct::pi * pow(max_r, 3.0);
  f1->SetParameter(0, area);

  std::vector<std::vector<double>> in{
      {1.4999, 0.0, 0.0},      {0.0, 1.49999, 0.0}, {0.0, 0.0, 1.4999},
      {1.4999, 0.000001, 0.0}, {1.0, 1.0, 0.0},     {0.74, 0.74, 0.74}};
  std::vector<std::vector<double>> out{{1.50001, 0.0, 0.0},
                                       {0.0, 1.50001, 0.0},
                                       {0.0, 0.0, 1.50001},
                                       {1.4999, 0.5, 0.0}};

  for (auto &e : in) {
    EXPECT_NEAR(f1->Eval(e[0], e[1], e[2]), 1.0 / volume, 1e-5);
  }
  for (auto &e : out) {
    EXPECT_NEAR(f1->Eval(e[0], e[1], e[2]), 0.0, 1e-5);
  }
}

TEST(function, step_function_3d_params) {
  EXPECT_EQ(sct::StepFunction_npar, 1);
  std::pair<unsigned, double> true_val = {0, 1.0};
  std::pair<unsigned, double> lookup_val = sct::StepFunction_params.at("sigma");
  EXPECT_EQ(true_val, lookup_val);
}

TEST(function, gaussian_params) {
  EXPECT_EQ(sct::Gaussian_npar, 1);
  std::pair<unsigned, double> true_val = {0, 1.0};
  std::pair<unsigned, double> lookup_val = sct::Gaussian_params.at("sigma");
  EXPECT_EQ(true_val, lookup_val);
}

TEST(function, WS_spherical_params) {
  EXPECT_EQ(sct::WoodsSaxonSpherical_npar, 2);
  std::pair<unsigned, double> true_val_r = {0, 5.0};
  std::pair<unsigned, double> true_val_skin_depth = {1, 0.5};
  std::pair<unsigned, double> lookup_val_r =
      sct::WoodsSaxonSpherical_params.at("radius");
  std::pair<unsigned, double> lookup_val_skin_depth =
      sct::WoodsSaxonSpherical_params.at("skin_depth");
  EXPECT_EQ(true_val_r, lookup_val_r);
  EXPECT_EQ(true_val_skin_depth, lookup_val_skin_depth);
}

TEST(function, WS_deformed_params) {
  EXPECT_EQ(sct::WoodsSaxonDeformed_npar, 4);
  std::pair<unsigned, double> true_val_r = {0, 5.0};
  std::pair<unsigned, double> true_val_skin_depth = {1, 0.5};
  std::pair<unsigned, double> true_val_b2 = {2, 0.1};
  std::pair<unsigned, double> true_val_b4 = {3, 0.1};
  std::pair<unsigned, double> lookup_val_r =
      sct::WoodsSaxonDeformed_params.at("radius");
  std::pair<unsigned, double> lookup_val_skin_depth =
      sct::WoodsSaxonDeformed_params.at("skin_depth");
  std::pair<unsigned, double> lookup_val_b2 =
      sct::WoodsSaxonDeformed_params.at("beta2");
  std::pair<unsigned, double> lookup_val_b4 =
      sct::WoodsSaxonDeformed_params.at("beta4");
  EXPECT_EQ(true_val_r, lookup_val_r);
  EXPECT_EQ(true_val_skin_depth, lookup_val_skin_depth);
  EXPECT_EQ(true_val_b2, lookup_val_b2);
  EXPECT_EQ(true_val_b4, lookup_val_b4);
}

TEST(function, hulthen_params) {
  EXPECT_EQ(sct::HulthenPDF_npar, 2);
  std::pair<unsigned, double> true_val_a = {0, 0.228};
  std::pair<unsigned, double> true_val_b = {1, 1.177};
  std::pair<unsigned, double> lookup_val_a = sct::HulthenPDF_params.at("a");
  std::pair<unsigned, double> lookup_val_b = sct::HulthenPDF_params.at("b");
  EXPECT_EQ(true_val_a, lookup_val_a);
  EXPECT_EQ(true_val_b, lookup_val_b);
}