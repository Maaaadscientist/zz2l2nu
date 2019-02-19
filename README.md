Instruction to produce the plots:
=================================

**a) The main analysis**

Everything can be launched from the launchAnalysis script. To know the options of this script, just run:
```
sh launchAnalysis.sh
```

If you want to *choose where the jobs will be saved (and avoir erasing your previous results)*, open this file and change the `suffix` of the analysis you want to run.
If you want to *change the jobs that will be ran on*, have a look at the `listDataset` and the corresponding files. Be careful that there are name conventions. Try to mimick as much as possible what you already see for the names.

**b) The InstrMET building**

Weights for the InstrMET can be found in `WeightsAndDatadriven/InstrMET/`. If you want to reproduce them, just do:
```
sh computeInstrMETWeights.sh
```
If you want to perform the closure test, launch:
```
sh doClosureTest.sh
```

**c) Running systematics**

A complete list of all the systematics can be found in the file `systList.txt`. The first column is the precise name of the systematic, and next columns are wildcards to designate on which samples these systematics need to be applied.3 different modes can be used to run on the systematics (the default option is not to compute them):

1) Run on only 1 systematic, up or down, for example `ewk_up`:
```
sh launchAnalysis.sh --syst ewk_up 1
```
Be careful to respect the exact name of the systematic (it must appear in one of the lines of the `systList.txt` file). This will launch the analysis only with this systematic, without the nominal shape, and only on catalogs containing the content of at least one of the columns of the `systList.txt` file for this syst in their name (with this example, `WZTo3LNu` or `ZZTo2L2Nu`.

2) Run on 1 systematic and produce both the up and down shapes, and the nominal (for example with `ewk_up`):
```
sh launchAnalysis.sh --syst ewk 1
```
This will run on the nominal shape and on `ewk_up` and `ewk_down`. If you don't want the nominal shapes to be produced, either comment the relevant line in the python file, or launch the script separately on the up and down shapes.

3) Run on all systematics and on nominal:
```
sh launchAnalysis.sh --syst all 1
```
This will run on all the systematics defined in the `systList.txt` file (if you want only some of these, you may comment some lines of this file). Caution: this command is meant to produce all the final shapes (to give to the datacards). Hence, by default, only the final results (mT shapes) and some necessary plots are kept, while all the control plots are dropped to make the output more readable. If you don't want that, change the relevant line in the python file.

Notice that all the plots produced with a given syst will have this syst in their name.

Later steps (2, 3 and 4) work in the same way as the standard analysis. In particular, notice that step 2 (harvesting) will also produce `_final` ROOT files for each type of samples; these files contain all the systs that were run in this production. Moreover, if you run with option `all`, the default behavior is that only these files are kept (without the `_final` in their name), in order to reduce the size of the outputs.

Just tell me how to runs the latest version of the code!
=========================================================

**Example 1: run on datadriven**
```
sh launchAnalysis.sh 0 #cleaning
sh launchAnalysis.sh 1 #run jobs
sh launchAnalysis.sh 2 #harvest jobs
sh launchAnalysis.sh 3 #run dataMCcomparison
sh launchAnalysis.sh 4 #publish results
```

**Example 2: run on MC-based analysis**
```
sh launchAnalysis.sh HZZanalysis 0 #cleaning
sh launchAnalysis.sh HZZanalysis 1 #run jobs
sh launchAnalysis.sh HZZanalysis 2 #harvest jobs
sh launchAnalysis.sh HZZanalysis 3 #run dataMCcomparison
sh launchAnalysis.sh HZZanalysis 4 #publish results
```

**Example 3: run on a specific bonzai**
   - create your datasetList and reference it in `launchAnalysis.sh` at the `listDataset` variable of the analysis you want to run.
   - check the naming convention! The better would probably to use the exact same name of what is found in the current list. (so if your sample is named "myDYsample" rename it "Bonzais-DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-001-ZZ2l2vPruner-MC_DLep_newTrigger.txt")
   - then launch the steps above

**Example 4: run on 10000 events**
   - open `runHZZ2l2nu.cc` and change the default value for `maxEvents`
   - then  launch the steps above

**Example 5: run on all systematics**
```
sh launchAnalysis.sh --syst all 0 #cleaning
sh launchAnalysis.sh --syst all 1 #run jobs
sh launchAnalysis.sh --syst all 2 #harvest jobs
```

Steps 3 and 4 can be run also, but will not do anything with the systematics (only the nominal).

**Example 6: run on a specific systematic**
```
sh launchAnalysis.sh --syst MYSYST 0 #cleaning
sh launchAnalysis.sh --syst MYSYST 1 #run jobs
sh launchAnalysis.sh --syst MYSYST 2 #harvest jobs
```

where `MYSYST` is one of the systs listed in the `systList.txt` file, without `_up` or `_down`.

Interacting with the IIHE cluster:
==============================
**1) Check my jobs**
```qstat -u $USER```

**2) Delete jobs**
```for j in $(qselect -u $USER);do timeout 3 qdel -a $j;done```

**3) Find jobs name and info**
```qstat -f $(qselect -u $USER -s EHQRTW) | grep "Job Id\|Job_Name\|resources_used.walltime"```

**4) Go interactively on the express queue for some tests**
```qsub -I -q express```

**5) Peek at your running jobs**
   - connect to m machine with option `-A`
   - do `qpeek jobID`

How to add a new looper/analysis:
=================================
To be written once the code is in a more stable form...
