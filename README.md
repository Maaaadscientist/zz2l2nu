Instruction to produce the plots:
=================================

**a) launch the jobs**
```
./prepareAllJobs.py --listDataset listSamplesToRun.txt --suffix aTest
big-submission sendJobs_aTest.cmd
```

**b) harvest the jobs**
(after the jobs are completed)
```
./prepareAllJobs.py --listDataset listSamplesToRun.txt --suffix aTest --harvest
```

**c) draw the plots**
With the `dataMCcomparison.C` macro
```root -l -q -b dataMCcomparison.C("HZZanalysis", "aTest")```

Additional options:
===================
**1) launch InstrMET analysis**
On step a) (only) described above, just add `--doInstrMETAnalysis`
On step b), change "HZZanalysis" by "InstrMET"

**2) launch on the express queue**
If one has a small amount of jobs to run and if those runs are very fast, you may consider running those on the express queue.
To do so, just add `--express` to the command described on step a).
The queue express has no queue time but only accept to run 8 jobs in parrallel with a wallTime of 32min. A maximum of 100 jobs can be submitted at the same time per user.
If we happen to use this a lot, the IT's are willing to give us more jobs.

**3) cleaning the directory**
To clean all the inputs AND outputs of your submission, launch: `sh launchAnalysis.sh 0`.
Be very careful since this will ERASE EVERYTHING you produced in this folder!

**4) run on a single root file**
To be written and implemented...

Using the submission scripts:
=============================
To be written...


How to add a new looper/analysis:
=================================
To be written...
