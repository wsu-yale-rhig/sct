# STAR Centrality Tools (SCT)
### Nick Elsey (WSU)

Implementation of the STAR experiment's Monte-Carlo Glauber model  and associated centrality utilities.

## Modules

### Monte-Carlo Glauber

A Monte-Carlo Glauber generator simulates ion collisions by constructing nuclei using pre-defined PDFs. These randomly generated nuclei are then "overlayed" on top of each other with a given impact parameter, and collisions between individual nucleons of the two nuclei are counted. The sct::MCGlauber generator can produce collisions in symmetric or asymmetric systems, generating nuclei as small as a single proton or as large as Uranium (or even further, using user-defined inputs). Default parameters for nucleus shapes are provided for multiple species, including Au, Pb, U, Cu, d and p. Providing non-default parameters for these or other species is also relatively easy. The generator also has various utility options, such as running systematic variations of the input nucleus shape.

### Centrality Tools
Provides methods to model reference multiplicity from glauber distributions, and to calculate centrality definitions from the results. STAR uses a two-part multiplicity model that states that some fraction of the multiplicity measued in a collision comes from the number of participant nucleons of a collision (the soft production), and some comes from the number of binary collisions (hard production). The relative abundance of each is controlled by a parameters x (0 < x < 1). This two component multiplicity gives a number of "ancestor" particles, each of which then produces N final state particles, where N is sampled from a negative binomial distribution. Both x and the parameters of the negative binomial are parameters that are generally fit using optimization routines, but they are generally restricted due to physical models and prior measurements to exist in a given range. The multiplicity model also takes into account STAR's multiplicity dependent efficiency, as well as allows for the simulation of a trigger bias. The NBDFit class provides an optimization routine that will fit a simulated glauber refmult distribution to a given data refmult distribution, using a chi2 test. Once an optimal fit is found, tools for calculating the centrality bins are provided.

### Examples
Examples of how to use the various tools are available in sct/analysis. While these are called examples, they have been used to produce results that are currently being used in STAR.

## Dependencies

### required
ROOT (source and binaries available at https://root.cern.ch, as well as on github: https://github.com/root-project/root). If you want to build from source and don't have a 16 core workstation, be ready to take a nap :)

gflags and glog: both are available on multiple platform package managers, including homebrew and macports on mac. If you want to build from source, you can do:
```
git clone https://github.com/gflags/gflags.git && mkdir gflags/build && cd gflags/build && cmake -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=/path/to/install/loc .. && make install
git clone https://github.com/google/glog.git && mkdir glog/build && cd glog/build && cmake -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=/path/to/install/loc .. && make install
```
(Order matters - glog depends on gflags)


### optional
boost (for example analysis routines). If you want to build the test suite, sct will build the gtest and benchmark libraries internally and statically link the executables, so make sure to use the ``--recurse-submodules`` option when cloning.

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
