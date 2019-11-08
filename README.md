# H&rarr;ZZ&rarr;2&ell;2&nu; analysis

Doxygen documentation for C++ code is available [here](http://homepage.iihe.ac.be/~aapopov/hzz2l2nu/doc/). At present it is generated manually by Andrey and might not correspond to the latest version of the code.


## Computing environment and building

At the start of each session, set up the environment with

```sh
. ./env.sh
```

This script also stores the path to the base directory in environment variable `HZZ2L2NU_BASE`, which should then be used in scripts and compiled code to resolve relative paths to auxiliary data files, such as data-driven weights.

Build the package with the following commands:

```sh
mkdir build
cd build
cmake ..
make -j 4
```

The warning from CMake about the new version of Boost can be safely ignored. Executable `runHZZanalysis` is put into `$HZZ2L2NU_BASE/bin`, and it is accessible from `$PATH`. To rebuild the package after a change has been introduced to the code, repeat `make`. To start the build from scratch, remove the directory `build` and repeat the commands above.


## Running interactively

Computationally heavy part of the analysis is carried out by program `runHZZanalysis`. Here is an example command to run it interactively:

```sh
runHZZanalysis --config 2016.yaml --ddf /pnfs/iihe/cms/store/group/HZZ2l2nu/Production/2019-08-16_2016/DDF/ZZTo2L2Nu.yaml \
  --analysis Main --max-events 10000
```

The first parameter is the path to the master configuration file, such as [`2016.yaml`](config/2016.yaml). It provides global settings that affect all analyises and all datasets. The path is resolved with the help of [`FileInPath`](http://homepage.iihe.ac.be/%7Eaapopov/hzz2l2nu/doc/classFileInPath.html) service. Standard configuration files are located in directory `$HZZ2L2NU_BASE/config`, which is checked by `FileInPath` automatically.

The second parameter is the path to a [dataset definition file](https://gitlab.cern.ch/HZZ-IIHE/hzz2l2nu/wikis/dataset-definitions) (either a full one or a derived fragment). It provides paths to input files included in the dataset and all dataset-specific configuration parameters.

The last two parameters specify which analysis should be executed and the maximal number of events to process. A number of other command line parameters are supported, many of them also have shortcuts. The complete list can be obtained by running

```sh
runHZZanalysis --analysis <analysis> --help
```

Note that different analyses support different sets of parameters (hence the flag `--analysis` above; without it the help for the default analysis is printed).


## Full analysis chain

When running the full analysis chain over a complete collection of datasets, the program `runHZZanalysis` is executed by steering scripts under the hood.


### a) The main analysis

Everything can be launched from the `launchAnalysis.sh` script. To know the options of this script, just run:

```sh
./launchAnalysis.sh
```

When called with appropriate arguments, this script will create a default task directory for auxiliary and output files. The name for this directory can be specified with option `--task-dir`.


### b) The InstrMET building

Weights for the InstrMET can be found in `WeightsAndDatadriven/InstrMET/`. If you want to reproduce them, just do:

```sh
sh computeInstrMETWeights.sh
```

If you want to perform the closure test, launch:

```sh
sh doClosureTest.sh
```

### c) Running systematics

A complete list of all supported systematic uncertainties can be found in file [`config/syst.yaml`](config/syst.yaml). For each uncertainty it provides a sequence of masks that defines which datasets are affected by this uncertainty. The masks are checked against the names of the datasets, as specified in [dataset definition files](https://gitlab.cern.ch/HZZ-IIHE/hzz2l2nu/wikis/dataset-definitions). Three different modes can be used to run on the systematics (the default option is not to compute them):

1. Evaluate one systematic variation, for example `ewk_up`:

   ```sh
   ./launchAnalysis.sh --syst ewk_up 1
   ```

   Labels of systematic variation are obtained by adding a postfix `_up` or `_down` to the names of uncertainties in `config/syst.yaml`. This will launch the analysis only with the requested variation, without the nominal shape, and only on datasets whose names are matched to the corresponding masks (with this example, `WZTo3LNu` or `ZZTo2L2Nu`).

2. Evaluate both up and down variations for a single systematic uncertainty (for example with `ewk`):

   ```sh
   ./launchAnalysis.sh --syst ewk 1
   ```

3. Run on the nominal configuration and all systematic variations:

   ```sh
   ./launchAnalysis.sh --syst all 1
   ```

   This will evaluate all uncertainties defined in file `config/syst.yaml` (if you want only some of these, you may comment some lines of this file). Caution: this command is meant to produce all the final shapes (to give to the datacards). Hence, by default, only the final results (m<sub>T</sub> shapes) and some necessary plots are kept, while all the control plots are dropped to make the output more readable. If you don't want that, change the relevant line in the python file.

Notice that all the plots produced with a given `<syst>` will have this `<syst>` in their name.

Later steps (2, 3 and 4) work in the same way as the standard analysis. In particular, notice that step 2 (harvesting) will also produce `_final` ROOT files for each type of samples; these files contain all the systs that were run in this production. Moreover, if you run with option `all`, the default behavior is that only these files are kept (without the `_final` in their name), in order to reduce the size of the outputs.


## Examples

### Run data-driven analysis

```sh
./launchAnalysis.sh 0  # cleaning
./launchAnalysis.sh 1  # run jobs
./launchAnalysis.sh 2  # harvest jobs
./launchAnalysis.sh 3  # run dataMCcomparison
./launchAnalysis.sh 4  # publish results
```

### Run MC-based analysis

```sh
./launchAnalysis.sh HZZanalysis 0  # cleaning
./launchAnalysis.sh HZZanalysis 1  # run jobs
./launchAnalysis.sh HZZanalysis 2  # harvest jobs
./launchAnalysis.sh HZZanalysis 3  # run dataMCcomparison
./launchAnalysis.sh HZZanalysis 4  # publish results
```

### Run on a custom collection of datasets

* Create a text file with paths to desired dataset definition files, similar to files `listSamples*.txt`.
* Overwrite `$listDataset` variable(s) for the analysis you want to run in the configuration section in `launchAnalysis.sh`.
* Execute the commands from the examples above.

If the total number of events to be processed is a few millions or smaller, it might be easier to just run over the datasets in question interactively.

### Run on all systematics

```sh
./launchAnalysis.sh --syst all 0  # cleaning
./launchAnalysis.sh --syst all 1  # run jobs
./launchAnalysis.sh --syst all 2  # harvest jobs
```

Steps 3 and 4 can be run also, but will not do anything with the systematics (only the nominal).

### Run on a specific systematic

```sh
./launchAnalysis.sh --syst <syst> 0  # cleaning
./launchAnalysis.sh --syst <syst> 1  # run jobs
./launchAnalysis.sh --syst <syst> 2  # harvest jobs
```

where `<syst>` is one of the systs listed in the `systList.txt` file, without `_up` or `_down`.

### Run with MELA reweighted signal sample

```sh
./launchAnalysis.sh --syst all --mela 0  # cleaning
./launchAnalysis.sh --syst all --mela 1  # run jobs
./launchAnalysis.sh --syst all --mela 2  # harvest jobs
```

## Interacting with the IIHE cluster

1. Check my jobs

   ```sh
   qstat -u $USER
   ```

2. Delete jobs

   ```sh
   for j in $(qselect -u $USER); do timeout 3 qdel -a $j; done
   ```

3. Find jobs name and info

   ```sh
   qstat -f $(qselect -u $USER -s EHQRTW) | grep "Job Id\|Job_Name\|resources_used.walltime"
   ```

4. Go interactively on the express queue for some tests

   ```sh
   qsub -I -q express
   ```

5. Peek at your running jobs

   - Connect to an `m*` machine with option `-A`.
   - Do `qpeek <jobID>`.
