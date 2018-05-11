# STAR Centrality Tools (SCT)
### Nick Elsey (WSU)

Implementation of the STAR experiment's Monte-Carlo Glauber model  and associated centrality utilities.

## Structure

### Monte-Carlo Glauber
The MC generator can run either from pre-defined nuclei (Au197, for instance) or from user- specific inputs. If using pre-defined settings, the generator can also run systematic variations on the default parameters for systematic error analysis.

### Centrality Calculation
Provides methods to model reference multiplicity from glauber results. The goal of the centrality utilities is to generate a trigger and efficiency corrected reference multiplicity distribution, that can then be used to identify the centrality of a collision, since we can not directly measure the impact parameter. To do this, the multiplicity is modeled by randomly sampling the  NPart vs NColl distribution measured in the MC Glauber, and a two-component multiplicity model to generate an event initial multiplicity. For each "ancestor" particle then, a final multiplicity is sampled from a negative binomial. 
However, this means we have quite a few parameters to tune (x, the fraction of "hard" production in the two component model, the negative binomial parameters, and the efficiency parameters). These can be optimized using optimized a grid search over the parameter space, using a chi2 fit to the data refmult as the objective function.

## Dependencies

### required
ROOT (source and binaries available at https://root.cern.ch, as well as on github: https://github.com/root-project/root). If you want to build from source and don't have a 16 core workstation, be ready to take a nap :)

gflags and glog: both are available on multiple platform package managers, including homebrew and macports on mac. If you want to build from source, you can do:
```
git clone https://github.com/gflags/gflags.git && mkdir gflags/build && cd gflags/build && cmake .. && make install
git clone https://github.com/google/glog.git && mkdir glog/build && cd glog/build && cmake .. && make install
```
(Order matters - glog depends on gflags)


### optional
boost (for example analysis routines), gtest (for tests), and benchmark (for benchmarking routines)

## Build
Once ROOT, gflags and glog are installed on the user's system, building should be simple:
```
git clone --recurse-submodules https://github.com/nickelsey/sct.git
cd sct
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/location .. && make install
```
If you have the boost libraries, you can build the example analysis scripts with the cmake option
```
cmake -DBUILD_BINARIES=ON .. && make
```
If you want to run the tests and benchmark routines, the sct build script will build gtest and benchmark internally, so the user doesn't have to worry about it.
```
cmake -DBUILD_TEST=ON .. && make
```
