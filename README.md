coupling_xolotl
=====

This is a [MOOSE](https://mooseframework.inl.gov/getting_started/index.html) application wrapping [Xolotl](https://github.com/ORNL-Fusion/xolotl/wiki) a cluster dynamics code.

Here is how to install this application, if you would like to use a conda-based environment:

First setup a conda environment following the steps from [here](https://mooseframework.inl.gov/getting_started/installation/conda.html) up to "Install MOOSE". Instead, install a lower-level set of packages:

```bash
conda create -n coupling_xolotl
conda activate coupling_xolotl
conda install moose-mpi moose-tools boost=1.84.0
```

Then get the code:

```bash
git clone https://github.com/SciDAC-MOOSE-Xolotl-coupling-group/coupling_xolotl.git
cd coupling_xolotl
git submodule init
git submodule update
```

After obtaining the code and downloading the main dependencies, several configuration and build steps
need to be performed. In xolotl:

```bash
cd xolotl
mkdir build
cd build
cmake -DXolotl_BUILD_PETSC=ON -DXolotl_BUILD_HYPRE=ON -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=~/projects/coupling_xolotl/xolotl/install ../
make
make install
```

This assumes you have placed your `coupling_xolotl` app in a `~/projects` directory. Next, in MOOSE:

```bash
cd ../../moose
export PETSC_DIR=~/projects/coupling_xolotl/xolotl/build/external/petsc_install
MOOSE_JOBS=8 ./scripts/update_and_rebuild_libmesh.sh
./scripts/update_and_rebuild_wasp.sh
```

Finally, the `coupling_xolotl` application can be built:

```bash
cd ..
make
```

Tests can then be run by calling the test script:

```bash
./run_tests
```

If your machine has N cores available the installation can go faster by using:
```bash
make -j N
```

If you have your own Boost installation you can define `BOOST_ROOT` before starting the installation.

Troubleshouting
------

**Clang and OpenMP**

If you use Clang with OpenMP support but still get an error in the build libmesh step stating that your compiler does not support OpenMP, try setting:
```
export OPENMP_CXXFLAGS='-Xpreprocessor -fopenmp -lomp'
export FFLAGS='-L/usr/local/lib'
```
