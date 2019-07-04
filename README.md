# H&rarr;ZZ&rarr;2&ell;2&nu; analysis

Doxygen documentation for C++ code is available [here](http://homepage.iihe.ac.be/~aapopov/hzz2l2nu/doc/). At present it is generated manually by Andrey and might not correspond to the latest version of the code.


## Computing environment and building

Set up the environment with

```sh
. ./env.sh
```

CMSSW environment is not needed any more. This script also stores the path to the base directory in environment variable `HZZ2L2NU_BASE`, which should then be used in scripts and compiled code to resolve relative paths to auxiliary data files, such as data-driven weights.

Build the package with the following commands:

```sh
mkdir build
cd build
cmake ..
make -j 4
```

The warning from CMake about the new version of Boost can be safely ignored. Executable `runHZZanalysis` is put into `$HZZ2L2NU_BASE/bin`, and it is accessible from `$PATH`. To rebuild the package after a change has been introduced to the code, repeat `make`. To start the build from scratch, remove the directory `build` and repeat the commands.


## Instructions to produce the plots


### a) The main analysis

Everything can be launched from the `launchAnalysis.sh` script. To know the options of this script, just run:

```sh
./launchAnalysis.sh
```

If you want to *choose where the jobs will be saved (and avoid erasing your previous results)*, open this file and change the `suffix` of the analysis you want to run.
If you want to *change the files that will be ran on*, have a look at the `listDataset*` variables and the corresponding text files. Be careful that there are name conventions. Try to mimick as much as possible what you already see for the names.

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

A complete list of all the systematics can be found in the file `systList.txt`. The first column is the precise name of the systematic, and next columns are wildcards to designate on which samples these systematics need to be applied. Three different modes can be used to run on the systematics (the default option is not to compute them):

1. Run on only 1 systematic, up or down, for example `ewk_up`:

   ```sh
   ./launchAnalysis.sh --syst ewk_up 1
   ```

   Be careful to respect the exact name of the systematic (it must appear in one of the lines of the `systList.txt` file). This will launch the analysis only with this systematic, without the nominal shape, and only on catalogs containing the content of at least one of the columns of the `systList.txt` file for this syst in their name (with this example, `WZTo3LNu` or `ZZTo2L2Nu`.

2. Run on 1 systematic and produce both the up and down shapes, and the nominal (for example with `ewk_up`):

   ```sh
   ./launchAnalysis.sh --syst ewk 1
   ```

   This will run on the nominal shape and on `ewk_up` and `ewk_down`. If you don't want the nominal shapes to be produced, either comment the relevant line in the python file, or launch the script separately on the up and down shapes.

3. Run on all systematics and on nominal:

   ```sh
   ./launchAnalysis.sh --syst all 1
   ```

   This will run on all the systematics defined in the `systList.txt` file (if you want only some of these, you may comment some lines of this file). Caution: this command is meant to produce all the final shapes (to give to the datacards). Hence, by default, only the final results (mT shapes) and some necessary plots are kept, while all the control plots are dropped to make the output more readable. If you don't want that, change the relevant line in the python file.

Notice that all the plots produced with a given `<syst>` will have this `<syst>` in their name.

Later steps (2, 3 and 4) work in the same way as the standard analysis. In particular, notice that step 2 (harvesting) will also produce `_final` ROOT files for each type of samples; these files contain all the systs that were run in this production. Moreover, if you run with option `all`, the default behavior is that only these files are kept (without the `_final` in their name), in order to reduce the size of the outputs.


## Just tell me how to run the latest version of the code!

### Example 1: Run data-driven analysis

```sh
./launchAnalysis.sh 0  # cleaning
./launchAnalysis.sh 1  # run jobs
./launchAnalysis.sh 2  # harvest jobs
./launchAnalysis.sh 3  # run dataMCcomparison
./launchAnalysis.sh 4  # publish results
```

### Example 2: Run MC-based analysis

```sh
./launchAnalysis.sh HZZanalysis 0  # cleaning
./launchAnalysis.sh HZZanalysis 1  # run jobs
./launchAnalysis.sh HZZanalysis 2  # harvest jobs
./launchAnalysis.sh HZZanalysis 3  # run dataMCcomparison
./launchAnalysis.sh HZZanalysis 4  # publish results
```

### Example 3: Run on a specific bonzai

   - Create a text file with a list of files, similar to `listSamples*.txt`, and in `launchAnalysis.sh` assign a path to it to the `listDataset` variable of the analysis you want to run.
   - Check the naming convention! The better would probably to use the exact same name of what is found in the current list. (So if your sample is named "myDYsample" rename it "Bonzais-DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-001-ZZ2l2vPruner-MC_DLep_newTrigger.txt".)
   - Then launch the steps above

### Example 4: Run on 10000 events

   - open `runHZZ2l2nu.cc` and change the default value for `maxEvents`
   - then  launch the steps above

### Example 5: Run on all systematics

```sh
./launchAnalysis.sh --syst all 0  # cleaning
./launchAnalysis.sh --syst all 1  # run jobs
./launchAnalysis.sh --syst all 2  # harvest jobs
```

Steps 3 and 4 can be run also, but will not do anything with the systematics (only the nominal).

### Example 6: Run on a specific systematic

```sh
./launchAnalysis.sh --syst <syst> 0  # cleaning
./launchAnalysis.sh --syst <syst> 1  # run jobs
./launchAnalysis.sh --syst <syst> 2  # harvest jobs
```

where `<syst>` is one of the systs listed in the `systList.txt` file, without `_up` or `_down`.

### Example for MELA: Run with MELA reweighted signal sample

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


## How to add a new looper/analysis

To be written once the code is in a more stable form&hellip;
