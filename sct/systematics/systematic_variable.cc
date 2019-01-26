// sct/systematics/systematic_variable.cc

#include "sct/systematics/systematic_variable.hh"

#include "sct/core/logging.hh"
#include "sct/utils/histogram_info.hh"

#include "TDirectory.h"


 namespace sct {

    SystematicVariable::SystematicVariable() : y_(), name_(""), label_(""), bins_(0), low_(0), high_(0), cumulant_flag_(false) {}


    SystematicVariable::SystematicVariable(GlauberObservable obs) {
      init(obs);
      cumulant_flag_ = false;
    }
    

    void SystematicVariable::init(GlauberObservable obs) {
      clear();

      // Use a dictionary to lookup the binning & naming parameters
      HistogramInfo& histInfo = HistogramInfo::instance();
      y_ = obs;
      name_ = histInfo.name(obs);
      label_ = histInfo.label(obs);
      bins_ = histInfo.bins(obs);
      low_ = histInfo.lowEdge(obs);
      high_ = histInfo.highEdge(obs);
    }
    
    void SystematicVariable::add(GlauberMod  mod) {
      // first, check to make sure it hasn't been added already - 
      // SystematicVariable does not support multiple copies of the same modification
      
      if (modifications_.find(mod) != modifications_.end()) {
        LOG(ERROR) << "glauber modification: " << glauberModToString[mod] 
                   << " has previously been added to this Systematic variable," 
                   << " it has not been re-added or overwritten";
        return;
      }

      // add to list - so we have an easy way of knowing
      // what modifications were added.
      modifications_.insert(mod);

      // build base name
      string base_name =  glauberModToString[mod] + name_;

      // Use a dictionary to lookup the binning & naming parameters
      HistogramInfo& histInfo = HistogramInfo::instance();

      // build histograms
      for (auto& obs : x_axis_variables_) {
        // build unique names and titles. This assumes multiple SystematicVariables
        // across the same root namespace do not have the same name, or there will
        // be name collisions.
        string name = base_name + histInfo.name(obs);
        string name_1d = name + "1d";
        string name_2d = name + "2d";
        string name_prof = name + "prof";
        string name_weight = name + "weight";
        string title = MakeString(";", histInfo.label(obs), ";", name_);
        string title_weight = MakeString(";", histInfo.label(obs), ";weight");
        
        // get x axis binning parameters
        unsigned xbins = histInfo.bins(obs);
        double xlow = histInfo.lowEdge(obs);
        double xhigh = histInfo.highEdge(obs);

        th1_.add(name_1d, title, xbins, xlow, xhigh);
        tprof_.add(name_prof, title, xbins, xlow, xhigh);
        th2_.add(name_2d, title, xbins, xlow, xhigh, bins_, low_, high_);
        weights_.add(name_weight, title_weight, xbins, xlow, xhigh);

        // check if we are calculating cumulants for this observable
        if (cumulant_flag_) {
          for (auto& order : cumulant_order_) {
            string name_1d_moment = name_1d + "moment" + std::to_string(order);
            string name_prof_moment = name_prof + "moment" + std::to_string(order);
            string name_1d_cumulant = name_1d + "cumulant" + std::to_string(order);

            string title_cumulant = MakeString(";", histInfo.label(obs), ";", name_);

            moment_th1_ .add(name_1d_moment, title_cumulant, xbins, xlow, xhigh);
            moment_tprof_.add(name_prof_moment, title_cumulant, xbins, xlow, xhigh);

          }
        }
      }
    }

    // delete all histograms
    void SystematicVariable::clear() {
      modifications_.clear();
      th1_.clear();
      th2_.clear();
      tprof_.clear();
      weights_.clear();
      moment_th1_ .clear();
      moment_tprof_.clear();
      cumulant_th1_.clear();
    }
    
    bool SystematicVariable::fillEvent(GlauberMod mod, EventDict<double>& dict, double weight) {
      // check to make sure GlauberMod was added

      // loop over all x axis variables
      for (auto& x : x_axis_variables_) {
        double x_val = dict[x];
        double y_val = dict[y_];

        // build names for 2d, prof and weight histograms
        string name_2d = getHistogramName(mod, x, "2d");
        string name_prof = getHistogramName(mod, x, "prof");
        string name_weight = getHistogramName(mod, x, "weight");

        th2_.fill(name_2d, x_val, y_val);
        tprof_.fill(name_prof, x_val, y_val * weight);
        weights_.fill(name_weight, x_val, weight);

        // fill moments if asked for
        if (cumulant_flag_) {
          for (auto& order : cumulant_order_) {
            // get historam names
            string name_moment = name_prof + MakeString("moment", std::to_string(order));
            double weighted_moment = pow(y_val, order) * weight;
            moment_tprof_.fill(name_moment, x_val, weighted_moment);
          }
        }
      }
     return true;
    }

    void SystematicVariable::write(TFile* file) {
      
      // get current directory so we can cd() back
      TDirectory* current_dir = TDirectory::CurrentDirectory();

      // change current directory to file
      file->cd();

      // start to write out, filling in the TH1Ds with proper weighted results as we go
      for (auto& mod : modifications_) {
        for (auto& x : x_axis_variables_) {

          // TH1 contents should be the weighted entries, profile / weight
          // we reweight by hand to get this
          TProfile* tmp_p = tprof_.get(getHistogramName(mod, x, "prof"));
          TProfile* tmp_weight = weights_.get(getHistogramName(mod, x, "weight"));
          TH1D* hist_1d = th1_.get(getHistogramName(mod, x, "1d"));
          reweight(tmp_p, tmp_weight, hist_1d);

          
          if (cumulant_flag_ && cumulant_order_.size() > 0) {
            // if we have moment histograms, perform reweighting for these as well.
            // collect pointers to the corrected moments and the cumulantsfor easy
            //  access when calculating the cumulants
            std::vector<TH1D*> corrected_moments(*cumulant_order_.rbegin() + 1, nullptr);
            std::vector<TH1D*> cumulants(*cumulant_order_.rbegin() + 1, nullptr);
            for (auto& order : cumulant_order_) {
              TProfile* tmp_moment_p = tprof_.get(getHistogramName(mod, x, MakeString("prof", "moment", std::to_string(order))));
              TH1D* tmp_moment_1d = th1_.get(getHistogramName(mod, x, MakeString("1d", "moment", std::to_string(order))));
              TH1D* tmp_cumulant_1d = th1_.get(getHistogramName(mod, x, MakeString("1d", "cumulant", std::to_string(order))));
              
              reweight(tmp_moment_p, tmp_weight, tmp_moment_1d);
              
              corrected_moments[order] = tmp_moment_1d;
              cumulants[order] = tmp_cumulant_1d;
            }

            // now we calculate the cumulants bin-by-bin
            for (int i = 1; i <=  corrected_moments[*cumulant_order_.begin()]->GetNbinsX(); ++i) {
              // we make vectors of the moments/errors. we need even moments up to the maximum value
              // specified in cumulant_order_. This is a std::set so the maximum value is always stored
              // last, so we get it with rbegin() 
              vector<double> moments(*cumulant_order_.rbegin()+1, 0);
              vector<double> moment_errors(*cumulant_order_.rbegin()+1, 0);

              // now each cumulant only depends on moments of the same order or less, so we can collect
              // the moments and calculate the cumulants in the same loop
              for (auto& order : cumulant_order_) {
                moments[order] = corrected_moments[order]->GetBinContent(i);
                moment_errors[order] = corrected_moments[order]->GetBinError(i);

                double cumulant = nthOrderCumulant(moments, order);
                double cumulant_error = nthOrderCumulantError(moments, moment_errors, order);

                cumulants[order]->SetBinContent(i, cumulant);
                cumulants[order]->SetBinError(i, cumulant_error);
              }
            }
          }
        }
      }

      // write to file
      th1_.write();
      th2_.write();
      tprof_.write();
      weights_.write();
      moment_th1_.write();
      moment_tprof_.write();

      // cd back to the directory root was in before
      current_dir->cd();
    }

    string SystematicVariable::getHistogramName(GlauberMod mod, GlauberObservable x_axis, string tag, unsigned cumulant_order) {
      return glauberModToString[mod] + name_ 
             + HistogramInfo::instance().name(x_axis) + tag 
             + (cumulant_order > 0 ? std::to_string(cumulant_order) : "");
    }

    void SystematicVariable::reweight(TProfile* source, TProfile* weight, TH1D* target) {
      // TProfile Division is not what we want - we need to create projections
      // into TH1Ds first
      unique_ptr<TH1D> tmp_p(source->ProjectionX());
      unique_ptr<TH1D> tmp_w(weight->ProjectionX());
      tmp_p->Divide(tmp_w.get());
      for (int i = 1; i <= target->GetNbinsX(); ++i) {
        target->SetBinContent(i, tmp_p->GetBinContent(i));
        target->SetBinError(i, tmp_p->GetBinError(i));
      }
    }

    double SystematicVariable::nthOrderCumulant(vector<double> moments, unsigned order) {
      // return an unphysical number if anything breaks
      double failure = 999999.0;

      // check to make sure we have the correct number of moments
      if (moments.size() < order) {
        LOG(ERROR) << "moments vector: " << moments << " does not contain enough moments "
                   << "to be calculating the " << order << " cumulant. Must have at least "
                   << order << " moments.";
        return failure;
      }

      double mu_2 = 0.0;
      double mu_4 = 0.0;
      double mu_6 = 0.0;

      switch (order) {
      case 2:
        mu_2 = moments[2];
        return mu_2 > 0.0 ? mu_2 : failure;
      
      case 4:
        mu_2 = moments[2];
        mu_4 = moments[4];
        return pow(abs(2.0 * pow(mu_2, 2) - mu_4), 1.0/4.0);
      
      case 6:
        mu_2 = moments[2];
        mu_4 = moments[4];
        mu_6 = moments[6];
        return pow(abs(0.25*(mu_6 - 9.0 * mu_4 * mu_2 + 12.0 * pow(mu_2, 3))), 1.0 / 6.0);
      
      default:
        LOG(ERROR) << "cumulant not implemented for order " << order;
        return failure;
      }
    }

    double SystematicVariable::nthOrderCumulantError(vector<double> moments, vector<double> errors, unsigned order) {

      // return an unphysical number if anything breaks
      double failure = 999999.0;

      // check to make sure we have the correct number of moments
      if (moments.size() < order) {
        LOG(ERROR) << "moments vector: " << moments << " does not contain enough moments "
                   << "to be calculating the " << order << " cumulant. Must have at least "
                   << order << " moments.";
        return failure;
      }

      // check to make sure we have the correct number of errors
      if (moments.size() < order) {
        LOG(ERROR) << "moments vector: " << errors << " does not contain enough errors "
                   << "to be calculating the " << order << " cumulant. Must have at least "
                   << order << " errors.";
        return failure;
      }

      // weighting defined by the order
      double order_weight = 1.0 / order;

      double mu_2 = 0.0, mu_4 = 0.0, mu_6 = 0.0;
      double mu_2_err = 0.0, mu_4_err = 0.0, mu_6_err = 0.0;
      double cumulant_4 = 0.0, error_4 = 0.0;
      double error_6_1 = 0.0, error_6_2 = 0.0, error_6 = 0.0, cumulant_6 = 0.0;
      
      switch (order) {
      case 2:
        mu_2 = moments[2];
        mu_2_err = errors[2];

        return mu_2_err * order_weight * pow(abs(mu_2), order_weight - 1.0);

      case 4:
        mu_2 = moments[2];
        mu_4 = moments[4];

        mu_2_err = errors[2];
        mu_4_err = errors[4];

        cumulant_4 = nthOrderCumulant(moments, 4);
        error_4 = sqrt(16.0 * pow(mu_2 * mu_2_err, 2.0) + mu_4_err * mu_4_err);
        return error_4 * order_weight * pow(abs(cumulant_4), order_weight - 1.0);
      
      case 6:
        mu_2 = moments[2];
        mu_4 = moments[4];
        mu_6 = moments[6];

        mu_2_err = errors[2];
        mu_4_err = errors[4];
        mu_6_err = errors[6];

        error_6_1  = mu_4 == 0.0 || mu_2 == 0.0 ? 0.0 :
                     9.0 * abs(mu_4*mu_2) * sqrt(pow(mu_4_err/mu_4, 2.0) + pow(mu_2_err/mu_2, 2.0));
        error_6_2  = 12.0 * 3.0 * mu_2 * mu_2 * mu_2_err;

        cumulant_6 = nthOrderCumulant(moments, 6);        
        error_6    = 0.25 * sqrt(pow(mu_6_err, 2.0) + pow(error_6_1, 2.0) + pow(error_6_2, 2.0));
        
        return error_6 * order_weight * pow(abs(cumulant_6), order_weight - 1.0);
      
      default:
        LOG(ERROR) << "cumulant not implemented for order " << order;
        return failure;
      }
    }

 } // namespace sct