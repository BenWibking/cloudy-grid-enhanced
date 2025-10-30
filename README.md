# Modified version of Cloudy

This version of Cloudy has been modified with various enhancements to produce
cooling tables for simulations under the assumption of ionization equilibrium.

## Example: isrf_ism cooling grid workflow

1. Build the solver if `cloudy.exe` is not present:

    ```bash
    make -C source -j4
    ```

2. Run the grid calculation from the repository root. This reads the native Cloudy
   input and writes `grid*_isrf_ism_*.txt` plus `isrf_ism_grid.grd`/`isrf_ism_summary.txt`
   in-place:

    ```bash
    ./source/cloudy.exe < isrf_ism_cooling_grid.in
    ```

3. Post-process the raw grid outputs into CIAOLoop-style ASCII tables and a stats
   summary (files land in `converted_tables/` by default):

    ```bash
    python3 scripts/postprocess_cooling_grid.py --prefix isrf_ism --output-dir converted_tables
    ```

4. Convert the ASCII tables into a single HDF5 cooling table. The converter infers
   the thermodynamic grid from the `.in` file; point it at the directory created in
   the previous step:

    ```bash
    python3 scripts/convert_cooling_grid_to_hdf5.py \
        isrf_ism_cooling_grid.in \
        converted_tables/isrf_ism.h5 \
        --data-dir converted_tables \
        --prefix isrf_ism
    ```

The resulting `converted_tables/isrf_ism.h5` collects the temperature, heating,
cooling, and mean molecular weight arrays in big-endian double precision. You
can inspect it with `h5dump` or load it directly in Python via `h5py`.

The original Cloudy README follows below.

# Cloudy

Cloudy is an _ab initio_ spectral synthesis code designed to model a wide range
of interstellar "clouds", from H II regions and planetary nebulae, to Active
Galactic Nuclei, and the hot intracluster medium that permeates galaxy clusters.

Cloudy has been in continuous development since 1978, led by Gary Ferland, and
in close collaboration with a number of scientists -- see the
[list of contributors](others.txt).


# Version

The current version of Cloudy is C25, released in 2025.
A summary of what is new is available
[here](https://gitlab.nublado.org/cloudy/cloudy/-/wikis/NewC25).

If you used Cloudy in your research, please cite our most recent
[release paper](https://ui.adsabs.harvard.edu/abs/2025arXiv250801102G)

## Brief History

Cloudy recently migrated to a pure git version control system from a
subversion (SVN) system (with limited support for git).
Cloudy had been on a SVN repository for about 15 years, which is still 
maintained as a read-only reference at
[https://trac.nublado.org](https://trac.nublado.org).

The migration was done on 2020 Dec 2 at r14364.
Only the trunk and a few actively maintained branches were migrated.

Previous releases of Cloudy are still available on the SVN site,
as well as tarballs on our
[release folder](https://data.nublado.org/cloudy_releases).


# Directory structure

There are seven directories below this one containing the:
1. atomic, molecular, grains data, as well as SEDs (```data/```);
1. documentation (```docs/```);
1. doxygen setup files (```doxygen/```);
1. a unit test library (```library/```);
1. some helpful scripts (```scripts/```);
1. the source (```source/```);
1. and test suite (```tsuite/```).
The test suite directory, tsuite, has a number of directories below it,
each exercising different aspects of Cloudy.

These directories contain all files needed to build and execute Cloudy.
Each directory has a readme file giving more information on its contents.
It is important to maintain this directory structure when the download is
opened on your computer.


# Building Cloudy

Instructions for building the code on various platforms are available on
[the wiki](https://gitlab.nublado.org/cloudy/cloudy/-/wikis/CompileCode).
Makefiles are provided for most popular compilers (see ```source/sys_*```).


# Documentation

The ```docs``` directory contains Hazy, Cloudy's documentation.


# API

Cloudy's API is described with [Doxygen](https://doxygen.nl).
A precompiled version is available
[online](https://data.nublado.org/doxygen/c22.00).

If you wish to produce a local instance, please follow the instructions in
the ```doxygen/``` directory.


# Contact us

See the project's [website](https://nublado.org) for new versions, bug fixes,
etc.

If you have any questions, please post them on the
[Cloudy user group](https://cloudyastrophysics.groups.io/g/Main/topics).
