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

TO BE WRITTEN BY NICOLAS

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
   - create yout datasetList and reference it in `launchAnalysis.sh` at the `listDataset` variable of the analysis you want to run.
   - check the naming convention! The better would probably to use the exact same name of what is found in the current list. (so if your sample is named "myDYsample" rename it "Bonzais-DYJetsToLL_M-50_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8-001-ZZ2l2vPruner-MC_DLep_newTrigger.txt")
   - then launch the steps above

**Example 4: run on 10000 events**
   - open `runHZZ2l2nu.cc` and change the default value for `maxEvents`
   - then  launch the steps above

**Example 5: run on all systematics**
TO BE WRITTEN BY NICOLAS

**Example 6: run on a specific systematic**
TO BE WRITTEN BY NICOLAS

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
