#ifndef SCT_SYSTEMATICS_GLAUBER_SYSTEMATICS_H
#define SCT_SYSTEMATICS_GLAUBER_SYSTEMATICS_H

#include "sct/centrality/multiplicity_model.h"
#include "sct/systematics/systematic_variable.h"
#include "sct/glauber/glauber_tree.h"
#include "sct/lib/memory.h"
#include "sct/lib/map.h"

#include <array>
#include <set>
#include <vector>

namespace sct {

  class GlauberSystematics {
  public:
    
    GlauberSystematics(const string& out_filename = "./");
    
    ~GlauberSystematics();

    bool setOutputFile(const string& out_file);
    
    // specify the input files that contain the glauber trees. If it has the 
    // systematic variation in the filename and it is one recognized by the sct
    // library, you don't need to manually provide the prefix. Otherwise, use
    // Load(filename, tag)
    bool load(const string& glauber_filename);
    bool load(const std::set<string>& glauber_filenames);
    bool load(const string& glauber_filename, GlauberMod tag);
    
    // need to set the multiplicity model parameters to the best-fit to data 
    // to generate a multiplicity on-the-fly for the glauber events
    void setMultiplicityModel(double npp, double k, double x, double ppEff,
                              double aaEff, double aaCent, double trigEff = 1.0,
                              bool constEff = false);

    // must supply centrality definition if one wants centrality results.
    // these should be lower bounds - i.e. 0-5%, multiplicity >= 520, 
    // 5-10% mult >= 480, etc. Requires a centrality definition that matches
    // the one defined in HistogramInfo. For instance, if HistogramInfo defines
    // a 17-bin definition, then cent_def should have 16 entries (one less,
    // because the last bin (usually 80-100%) is always defined from 0)
    void setCentralityBins(std::vector<double> cent_def);
    
    // flag to turn on/off npp variations in the multiplicity model as  another
    // variation. Looks for a variation named "nominal" when running, will error out
    // otherwise
    void doNppVariations(bool flag = true) {do_npp_variations_ = flag;}

    // use unit weight when calculating cumulants 
    void useUnitWeight(bool flag = true) {use_unit_weight_ = flag;}

    // runs systematic analysis over all inputs. Calls ProcessFile() for each 
    // valid TFile
    bool run();
    
    // writes results to disk in directory out_dir_
    bool write();

  private:
    
    // used to create all histograms and initialize output
    void initialize();

    // deletes all histograms
    void clearHistograms();

    // attempts to open a root file and add it to the dictionary
    bool loadFile(const string& filename, GlauberMod tag);

    // attempts to identify which modification was used when running
    // the glauber model from the file name and returns an ID as an enum
    GlauberMod findFileTag(const string& filename);

    // used to calculate centrality from the multiplicity
    int getCentrality(double multiplicity);

    // used to check that output file can be created
    bool openOutput();

    // runs a single file (+- npp variation)
    void runFile(GlauberMod key, shared_ptr<TFile> file);
    
    // output file and output file name
    string out_filename_;
    unique_ptr<TFile> out_file_;
    
    // holder for all TFiles (should be one file per variation)
    sct_map<GlauberMod, shared_ptr<TFile>, EnumClassHash> variations_;
    
    sct_map<GlauberMod, unique_ptr<MultiplicityModel>, EnumClassHash> mult_model_;

    // centrality lower bounds
    std::vector<double> centrality_lower_bounds_;

    // flags 
    bool do_npp_variations_;
    bool use_unit_weight_;
    
    // histogram containers
    SystematicVariable impact_parameter_;           // impact parameter
    SystematicVariable n_part_;                     // number of participants
    SystematicVariable n_coll_;                     // number of binary collisions
    SystematicVariable multiplicity_;               // event multiplicity (uncorrected)
    SystematicVariable rp_area_;                    // reaction plane area
    SystematicVariable pp_area_;                    // participant plane area
    SystematicVariable rp_ecc_;                     // npart weighted second order reaction plane eccentricity
    SystematicVariable rp_ecc_mult_;                // multiplicity weighted second order reaction plane eccentricity
    std::array<SystematicVariable, 3> pp_ecc_;      // npart weighted [second,third,fourth] order participant plane eccentricity
    std::array<SystematicVariable, 3> pp_ecc_mult_; // multiplicity weight
    
  };
  
} // namespace sct

#endif // SCT_SYSTEMATICS_GLAUBER_SYSTEMATICS_H
