// ugly hack to test private member SystematicVariable -
// namely, the cumulant and cumulant error calculations
#define protected public
#define private public
#include "sct/systematics/systematic_variable.h"
#undef protected
#undef private

#include "sct/lib/enumerations.h"
#include "sct/lib/logging.h"
#include "sct/lib/string/string_utils.h"
#include "sct/utils/histogram_info.h"

#include "gtest/gtest.h"

#include <cmath>
#include <set>
#include <string>

using std::string;

// test histograms are initialized properly when
// multiple GlauberMods are added

TEST(SystematicVariable, initialization) {
  std::set<sct::GlauberMod> mods{sct::GlauberMod::Nominal,
                                 sct::GlauberMod::Large,
                                 sct::GlauberMod::SmallXSec};

  std::set<sct::GlauberObservable> x_axis_observables{
      sct::GlauberObservable::B, sct::GlauberObservable::Ncoll,
      sct::GlauberObservable::Npart, sct::GlauberObservable::Multiplicity};

  std::set<sct::HistType> built_types{sct::HistType::TH1, sct::HistType::TH2,
                                      sct::HistType::Prof,
                                      sct::HistType::Weight};

  std::set<sct::HistType> empty_types{sct::HistType::Cumulant,
                                      sct::HistType::MomentProf,
                                      sct::HistType::MomentTH1};

  sct::SystematicVariable var(sct::GlauberObservable::Centrality);

  for (auto& mod : mods) {
    var.add(mod);
  }

  for (auto& mod : mods) {
    for (auto& x_axis : x_axis_observables) {
      for (auto& type : built_types) {
        EXPECT_NE(var.get(mod, x_axis, type), nullptr);
        for (auto& order : {2, 4, 6})
          EXPECT_EQ(var.get(mod, x_axis, type, order), nullptr);
      }

      for (auto& type : empty_types) {
        for (auto& order : {2, 4, 6})
          EXPECT_EQ(var.get(mod, x_axis, type, order), nullptr);
      }
    }
  }
}

// test to make sure cumulant initialization happens
// properly after calculateCumulants() is called

TEST(SystematicVariable, cumulantInitialization) {
  std::set<sct::GlauberMod> mods{sct::GlauberMod::Nominal,
                                 sct::GlauberMod::Large,
                                 sct::GlauberMod::SmallXSec};

  std::set<sct::GlauberObservable> x_axis_observables{
      sct::GlauberObservable::B, sct::GlauberObservable::Ncoll,
      sct::GlauberObservable::Npart, sct::GlauberObservable::Multiplicity};

  std::set<sct::HistType> built_types{sct::HistType::TH1, sct::HistType::TH2,
                                      sct::HistType::Prof,
                                      sct::HistType::Weight};

  std::set<sct::HistType> empty_types{sct::HistType::Cumulant,
                                      sct::HistType::MomentProf,
                                      sct::HistType::MomentTH1};

  sct::SystematicVariable var(sct::GlauberObservable::Centrality);
  var.calculateCumulants();

  for (auto& mod : mods) {
    var.add(mod);
  }

  for (auto& mod : mods) {
    for (auto& x_axis : x_axis_observables) {
      for (auto& type : built_types) {
        EXPECT_NE(var.get(mod, x_axis, type), nullptr);
        for (auto& order : {2, 4, 6})
          EXPECT_EQ(var.get(mod, x_axis, type, order), nullptr);
      }

      for (auto& type : cumulant_types) {
        EXPECT_EQ(var.get(mod, x_axis, type), nullptr);
        for (auto& order : {2, 4, 6})
          EXPECT_NE(var.get(mod, x_axis, type, order), nullptr);
      }

      for (auto& type : empty_types) {
        EXPECT_EQ(var.get(mod, x_axis, type), nullptr);
      }
    }
  }
}

// test to make sure filling routines are filling histograms
// as expected

TEST(SystematicVariable, filling) {
  std::set<sct::GlauberMod> mods{sct::GlauberMod::Nominal,
                                 sct::GlauberMod::Large,
                                 sct::GlauberMod::SmallXSec};

  std::set<sct::GlauberObservable> x_axis_observables{
      sct::GlauberObservable::B, sct::GlauberObservable::Ncoll,
      sct::GlauberObservable::Npart, sct::GlauberObservable::Multiplicity};

  std::set<sct::HistType> built_types{sct::HistType::TH1, sct::HistType::TH2,
                                      sct::HistType::Prof,
                                      sct::HistType::Weight};

  std::set<sct::HistType> empty_types{sct::HistType::Cumulant,
                                      sct::HistType::MomentProf,
                                      sct::HistType::MomentTH1};

  sct::SystematicVariable var(sct::GlauberObservable::Centrality);
  var.calculateCumulants();

  for (auto& mod : mods) {
    var.add(mod);
  }

  sct::EventDict<double> event;
  event[sct::GlauberObservable::Centrality] = 5;
  event[sct::GlauberObservable::Multiplicity] = 100;
  event[sct::GlauberObservable::Npart] = 23;
  event[sct::GlauberObservable::Ncoll] = 32;
  event[sct::GlauberObservable::B] = 0.1;

  double nominalweight = 1.0;
  var.fillEvent(sct::GlauberMod::Nominal, event, nominalweight);
  double largeweight = 5.0;
  var.fillEvent(sct::GlauberMod::Large, event, largeweight);

  for (auto& x_axis : x_axis_observables) {
    for (auto& mod : mods) {
      int nbins_x = sct::HistogramInfo::instance().bins(x_axis);
      double low_x = sct::HistogramInfo::instance().lowEdge(x_axis);
      double high_x = sct::HistogramInfo::instance().highEdge(x_axis);
      double delta_x = (high_x - low_x) / nbins_x;
      double x_val = event[x_axis];
      int bin_x = 0;
      for (int i = 0; i < nbins_x; ++i) {
        if (x_val < low_x + (delta_x * (i + 1))) {
          bin_x = i + 1;
          break;
        }
      }

      int nbins_y = sct::HistogramInfo::instance().bins(
          sct::GlauberObservable::Centrality);
      double low_y = sct::HistogramInfo::instance().lowEdge(
          sct::GlauberObservable::Centrality);
      double high_y = sct::HistogramInfo::instance().highEdge(
          sct::GlauberObservable::Centrality);
      double delta_y = (high_y - low_y) / nbins_y;
      double y_val = event[sct::GlauberObservable::Centrality];
      int bin_y = 0;
      for (int i = 0; i < nbins_y; ++i) {
        if (y_val < low_y + (delta_y * (i + 1))) {
          bin_y = i + 1;
          break;
        }
      }

      if (mod == sct::GlauberMod::Nominal) {
        EXPECT_EQ(var.get(mod, x_axis, sct::HistType::TH1)->Integral(), 0.0);
        EXPECT_EQ(var.get(mod, x_axis, sct::HistType::Prof)->Integral(),
                  nominalweight * event[sct::GlauberObservable::Centrality]);
        EXPECT_EQ(var.get(mod, x_axis, sct::HistType::Weight)->Integral(),
                  nominalweight);
        EXPECT_EQ(var.get(mod, x_axis, sct::HistType::TH2)->Integral(), 1.0);

        EXPECT_EQ(
            var.get(mod, x_axis, sct::HistType::Prof)->GetBinContent(bin_x),
            nominalweight * event[sct::GlauberObservable::Centrality]);
        EXPECT_EQ(
            var.get(mod, x_axis, sct::HistType::Weight)->GetBinContent(bin_x),
            nominalweight);

        int bin_2d =
            var.get(mod, x_axis, sct::HistType::TH2)->GetBin(bin_x, bin_y);
        EXPECT_EQ(
            var.get(mod, x_axis, sct::HistType::TH2)->GetBinContent(bin_2d),
            1.0);

        // test moments
        for (auto& order : {2, 4, 6}) {
          EXPECT_EQ(
              var.get(mod, x_axis, sct::HistType::MomentTH1, order)->Integral(),
              0.0);
          EXPECT_EQ(var.get(mod, x_axis, sct::HistType::MomentProf, order)
                        ->Integral(),
                    nominalweight *
                        pow(event[sct::GlauberObservable::Centrality], order));
          EXPECT_EQ(var.get(mod, x_axis, sct::HistType::MomentProf, order)
                        ->GetBinContent(bin_x),
                    nominalweight *
                        pow(event[sct::GlauberObservable::Centrality], order));
          EXPECT_EQ(
              var.get(mod, x_axis, sct::HistType::Cumulant, order)->Integral(),
              0.0);
        }
      }

      else if (mod == sct::GlauberMod::Large) {
        EXPECT_EQ(var.get(mod, x_axis, sct::HistType::TH1)->Integral(), 0.0);
        EXPECT_EQ(var.get(mod, x_axis, sct::HistType::Prof)->Integral(),
                  largeweight * event[sct::GlauberObservable::Centrality]);
        EXPECT_EQ(var.get(mod, x_axis, sct::HistType::Weight)->Integral(),
                  largeweight);
        EXPECT_EQ(var.get(mod, x_axis, sct::HistType::TH2)->Integral(), 1.0);

        EXPECT_EQ(
            var.get(mod, x_axis, sct::HistType::Prof)->GetBinContent(bin_x),
            largeweight * event[sct::GlauberObservable::Centrality]);
        EXPECT_EQ(
            var.get(mod, x_axis, sct::HistType::Weight)->GetBinContent(bin_x),
            largeweight);

        int bin_2d =
            var.get(mod, x_axis, sct::HistType::TH2)->GetBin(bin_x, bin_y);
        EXPECT_EQ(
            var.get(mod, x_axis, sct::HistType::TH2)->GetBinContent(bin_2d),
            1.0);

        // test moments
        for (auto& order : {2, 4, 6}) {
          EXPECT_EQ(
              var.get(mod, x_axis, sct::HistType::MomentTH1, order)->Integral(),
              0.0);
          EXPECT_EQ(var.get(mod, x_axis, sct::HistType::MomentProf, order)
                        ->Integral(),
                    largeweight *
                        pow(event[sct::GlauberObservable::Centrality], order));
          EXPECT_EQ(var.get(mod, x_axis, sct::HistType::MomentProf, order)
                        ->GetBinContent(bin_x),
                    largeweight *
                        pow(event[sct::GlauberObservable::Centrality], order));
          EXPECT_EQ(
              var.get(mod, x_axis, sct::HistType::Cumulant, order)->Integral(),
              0.0);
        }
      }

      else {
        EXPECT_EQ(var.get(mod, x_axis, sct::HistType::TH1)->Integral(), 0.0);
        EXPECT_EQ(var.get(mod, x_axis, sct::HistType::Prof)->Integral(), 0.0);
        EXPECT_EQ(var.get(mod, x_axis, sct::HistType::Weight)->Integral(), 0.0);
        EXPECT_EQ(var.get(mod, x_axis, sct::HistType::TH2)->Integral(), 0.0);

        for (auto& order : {2, 4, 6}) {
          EXPECT_EQ(var.get(mod, x_axis, sct::HistType::MomentProf, order)
                        ->Integral(),
                    0.0);
          EXPECT_EQ(
              var.get(mod, x_axis, sct::HistType::MomentTH1, order)->Integral(),
              0.0);
          EXPECT_EQ(
              var.get(mod, x_axis, sct::HistType::Cumulant, order)->Integral(),
              0.0);
        }
      }
    }
  }
}

// test that the cumulant calculated in sct is equal to that of
// the previous STAR glauber model

// old formulae taken from StGLauberCumulantHistogramMaker

double Get4thOrderCumulant(double c2, double c4);
double Get6thOrderCumulant(double c2, double c4, double c6);
double GetNthOrderCumulantError(double order, double val, double err);
double Get2ndOrderCumulantError(double c2, double c2error);
double Get4thOrderCumulantError(double c2, double c4, double c2error,
                                double c4error);
double Get6thOrderCumulantError(double c2, double c4, double c6, double c2error,
                                double c4error, double c6error);
double GetCumulant(unsigned order, double* val);
double GetCumulantError(unsigned order, double* val, double* err);

TEST(SystematicVariable, 2ndOrderCumulant) {
  sct::SystematicVariable var(sct::GlauberObservable::Centrality);

  std::vector<double> moments{0.0, 0.0, 2, 0.0, 2.5, 0.0, 1.0};
  double moments_old[] = {2, 2.5, 1.0};
  std::vector<double> errors{0.0, 0.0, 0.92, 0.0, 5.0, 0.0, 0.2};
  double errors_old[] = {0.92, 5.0, 0.2};

  double c2 = var.nthOrderCumulant(moments, 2);
  double c2_old = GetCumulant(2, moments_old);

  double c2_err = var.nthOrderCumulantError(moments, errors, 2);
  double c2_err_old = GetCumulantError(2, moments_old, errors_old);

  EXPECT_NEAR(c2, c2_old, 1e-3);
  EXPECT_NEAR(c2_err, c2_err_old, 1e-3);
}

TEST(SystematicVariable, 4thOrderCumulant) {
  sct::SystematicVariable var(sct::GlauberObservable::Centrality);

  std::vector<double> moments{0.0, 0.0, 2, 0.0, 2.5, 0.0, 1.0};
  double moments_old[] = {2, 2.5, 1.0};
  std::vector<double> errors{0.0, 0.0, 0.92, 0.0, 5.0, 0.0, 0.2};
  double errors_old[] = {0.92, 5.0, 0.2};

  double c4 = var.nthOrderCumulant(moments, 4);
  double c4_old = GetCumulant(4, moments_old);

  double c4_err = var.nthOrderCumulantError(moments, errors, 4);
  double c4_err_old = GetCumulantError(4, moments_old, errors_old);

  EXPECT_NEAR(c4, c4_old, 1e-3);
  EXPECT_NEAR(c4_err, c4_err_old, 1e-3);
}

TEST(SystematicVariable, 6thOrderCumulant) {
  sct::SystematicVariable var(sct::GlauberObservable::Centrality);

  std::vector<double> moments{0.0, 0.0, 2, 0.0, 2.5, 0.0, 1.0};
  double moments_old[] = {2, 2.5, 1.0};
  std::vector<double> errors{0.0, 0.0, 0.92, 0.0, 5.0, 0.0, 0.2};
  double errors_old[] = {0.92, 5.0, 0.2};

  double c6 = var.nthOrderCumulant(moments, 6);
  double c6_old = GetCumulant(6, moments_old);

  double c6_err = var.nthOrderCumulantError(moments, errors, 6);
  double c6_err_old = GetCumulantError(6, moments_old, errors_old);

  EXPECT_NEAR(c6, c6_old, 1e-3);
  EXPECT_NEAR(c6_err, c6_err_old, 1e-3);
}

double Get4thOrderCumulant(double c2, double c4) {
  double raw = 2.0 * c2 * c2 - c4;
  double sign = (raw > 0.0) ? 1.0 : -1.0;

  return pow(sign * raw, 0.25);
}

//____________________________________________________________________________________________________
double Get6thOrderCumulant(double c2, double c4, double c6) {
  double raw = 0.25 * (c6 - 9.0 * c4 * c2 + 12.0 * c2 * c2 * c2);
  double sign = (raw > 0.0) ? 1.0 : -1.0;

  return pow(sign * raw, 1.0 / 6.0);
}

//____________________________________________________________________________________________________
double GetNthOrderCumulantError(double order, double val, double err) {
  if (val == 0.0) return 0.0;

  return abs(order) * pow(abs(val), order - 1.0) * err;
}

//____________________________________________________________________________________________________
double Get2ndOrderCumulantError(double c2, double c2error) {
  return GetNthOrderCumulantError(0.5, c2, c2error);
}

//____________________________________________________________________________________________________
double Get4thOrderCumulantError(double c2, double c4, double c2error,
                                double c4error) {
  double cumulant = Get4thOrderCumulant(c2, c4);
  double error = sqrt(16.0 * c2 * c2 * c2error * c2error + c4error * c4error);

  return GetNthOrderCumulantError(0.25, cumulant, error);
}

//____________________________________________________________________________________________________
double Get6thOrderCumulantError(double c2, double c4, double c6, double c2error,
                                double c4error, double c6error) {
  double cumulant = Get6thOrderCumulant(c2, c4, c6);
  double error1 =
      (c4 != 0.0 && c2 != 0.0)
          ? 9.0 * abs(c4 * c2) *
                sqrt(pow(c4error / c4, 2.0) + pow(c2error / c2, 2.0))
          : 0.0;
  double error2 = 12.0 * 3.0 * c2 * c2 * c2error;
  double error =
      0.25 * sqrt(c6error * c6error + error1 * error1 + error2 * error2);

  return GetNthOrderCumulantError(1.0 / 6.0, cumulant, error);
}

//____________________________________________________________________________________________________
double GetCumulant(unsigned order, double* val) {
  /// val array should contain val[] = {c2, c4, c6, ...}
  if (!val) return -9999.;

  double invalid = -9999.;

  switch (order) {
    case 2:  // 2nd order cumulant
      return (val[0] < 0.0) ? invalid : sqrt(val[0]);

    case 4:  // 4th order cumulant
      return Get4thOrderCumulant(val[0], val[1]);

    case 6:  // 6th order cumulant
      return Get6thOrderCumulant(val[0], val[1], val[2]);

    default:
      assert(0);
  }

  return -9999.;
}

//____________________________________________________________________________________________________
double GetCumulantError(unsigned order, double* val, double* err) {
  /// val and err array should contain val[] = {c2, c4, c6, ...}, err[] =
  /// {c2error, c4error, c6error, ...}
  if (!val) return -9999.;
  if (!err) return -9999.;

  switch (order) {
    case 2:  // 2nd order cumulant
      return Get2ndOrderCumulantError(val[0], err[0]);

    case 4:  // 4th order cumulant
      return Get4thOrderCumulantError(val[0], val[1], err[0], err[1]);

    case 6:  // 6th order cumulant
      return Get6thOrderCumulantError(val[0], val[1], val[2], err[0], err[1],
                                      err[2]);

    default:
      assert(0);
  }

  return -9999.;
}