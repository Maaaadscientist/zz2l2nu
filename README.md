# H&rarr;ZZ&rarr;2&ell;2&nu; analysis

This repository uses Git submodules and hence should be cloned with

```sh
git clone --recurse-submodules ssh://git@gitlab.cern.ch:7999/HZZ-IIHE/hzz2l2nu.git
```

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
make -j $(nproc)
```

The warning from CMake about the new version of Boost can be safely ignored. Executable `runHZZanalysis` is put into `$HZZ2L2NU_BASE/bin`, and it is accessible from `$PATH`. To rebuild the package after a change has been introduced to the code, repeat `make`. To start the build from scratch, remove the directory `build` and repeat the commands above.

It is also possible to create a program outside of the repository and link it against the shared library of the framework. See [here](https://gitlab.cern.ch/HZZ-IIHE/hzz2l2nu/-/wikis/shared-library) for documentation.


## Running interactively

Computationally heavy part of the analysis is carried out by program `runHZZanalysis`. Here is an example command to run it interactively:

```sh
runHZZanalysis --config 2016.yaml \
  --ddf /pnfs/iihe/cms/store/group/HZZ2l2nu/Production/2020-04-16_2016-NanoAODv6/DDF/Dilepton/ZZTo2L2Nu.yaml \
  --analysis DileptonTrees --max-events 10000
```

The first parameter is the path to the master configuration file, such as [`2016.yaml`](config/2016.yaml). It provides global settings that affect all analyises and all datasets. The path is resolved with the help of [`FileInPath`](http://homepage.iihe.ac.be/%7Eaapopov/hzz2l2nu/doc/classFileInPath.html) service. Standard configuration files are located in directory `$HZZ2L2NU_BASE/config`, which is checked by `FileInPath` automatically.

The second parameter is the path to a [dataset definition file](https://gitlab.cern.ch/HZZ-IIHE/hzz2l2nu/wikis/dataset-definitions) (either a full one or a derived fragment). It provides paths to input files included in the dataset and all dataset-specific configuration parameters.

The last two parameters specify which analysis should be executed and the maximal number of events to process. A number of other command line parameters are supported, many of them also have shortcuts. The complete list can be obtained by running

```sh
runHZZanalysis --analysis <analysis> --help
```

Note that different analyses support different sets of parameters (hence the flag `--analysis` above; without it the help for the default analysis is printed).


## Using batch system

This repository includes scripts to submit jobs to the IIHE batch system and to merge their output.


### Submission

The submission is done with commands like

```sh
prepare_jobs.py --task-dir batch --config 2016.yaml \
  -- $HZZ2L2NU_BASE/config/samples_NRB_2016.txt -a DileptonTrees
cd batch
big-submission send_jobs.sh
```

The first parameter is the name of the directory in which job scripts and output files will be stored. As with the interactive running, the second parameter is the master configuration file; it is also automatically forwarded to `runHZZanalysis`. The first positional argument is a list of datasets to be processed. Each job will process at maximum one dataset, and large datasets are automatically split among multiple jobs. All subsequent positional arguments are forwarded to the executable without a change. When some of them start with dashes, as in the example above, it is convenient to denote the start of positional arguments with `--` so that they are not interpreted as optional arguments for `prepare_jobs.py`.

By default, `prepare_jobs.py` will run on nodes of the batch system executable `runHZZanalysis`. But it is possible to supply a different executable with option `--prog`.


### Systematic variations

The full list of options supported by `prepare_jobs.py` is available through its help message. One of them is option `--syst` that requests systematic variations to be applied. All supported systematic uncertainties are listed in file [`config/syst.yaml`](config/syst.yaml). For each uncertainty it provides a sequence of masks that defines which datasets are affected by this uncertainty. The masks are checked against the names of the datasets, as specified in [dataset definition files](https://gitlab.cern.ch/HZZ-IIHE/hzz2l2nu/wikis/dataset-definitions). There are several ways to specify a systematic variation:

* `--syst jec_up` or `--syst jec_down` <br /> These are examples of fully specified variations. They will be propagated to the underlying executable without a change, but only datasets that are affected by those variations (as indicated in `config/syst.yaml`) will be processed.
* `--syst jec` <br /> Jobs for both up and down variations of given type will be created.
* `--syst weights` <br /> Will evaluate all weight-based variations. Only applicable to analyses that produce trees with individual events.
* `--syst all` <br /> Jobs for all possible systematic variations, as well as the nominal configuration, will be created.


### Harvesting

When all jobs have finished (which can be checked with `qstat -u $USER`), their outputs can be merged with

```sh
harvest.py --task-dir batch --config 2016.yaml \
  $HZZ2L2NU_BASE/config/samples_NRB_2016.txt
```

The will be placed in `batch/merged`. You will probably want to move that directory somewhere and delete the rest of directory `batch`.


## Post-processing with event-based analysis

Instructions below are for the event-based analyses, such as `DileptonTrees`, only. They assume that the harvesting has finished, and the merged files are available in directory `$tree_dir`.

### Plots with comparison between data and simulation

These plots can be produced by running

```sh
plot_data_sim.py ${HZZ2L2NU_BASE}/config/plot_data_sim.yaml --prefix ${tree_dir}/ --output data_sim
```

Event selection and variables to be plotted are described in the configuration file given as the first argument to the script. The [file](config/plot_data_sim.yaml) included in the repository is considered an example, which you would adjust to your needs. Names (or rather name fragments) of input ROOT files are specified in the configuraiton. Each name is _textually_ prepended with the prefix given with the corresponding command line option  (or, alternatively, in the configuration file), which allows to specify the directory containing the files but also a common starting name fragment. In addition to this prefix, multiple standard endings are tried for each file name. The directory in which produced figures will be stored is given by flag `--output`.


### Templates for statistical analysis

Templates for statistical analysis can be constructed with

```sh
build_templates.py $tree_dir --output templates.root
```

This script will produce m<sub>T</sub> histograms. It expects that all systematic variations are available in the input files, so the should have been produced with `--syst all` option. Running this script takes around half an hour, so you should use one of the `mlong` user interfaces.

To visualize systematic variations in the produced templates, run

```sh
plot_syst_variations.py templates.py --output fig
```


## Full analysis chain

Script `launchAnalysis.sh` can execute the scripts discussed above for you, although at the price of a lower flexibility.


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

Request systematic variations with

```sh
./launchAnalysis.sh --syst <syst> 1
```

The value of `<syst>` is propagated to `prepare_jobs.py` discussed above. Notice that all the plots produced with a given `<syst>` will have this `<syst>` in their name.

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
