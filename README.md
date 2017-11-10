Instruction to produce the plots:
=================================

**a) launch the jobs**
```
./prepareAllJobs.py --listDataset listSamplesToRun.txt --suffix aTest
big-submission sendJobs_aTest
```

**b) harvest the jobs**
(after the jobs are completed)
```
./prepareAllJobs.py --listDataset listSamplesToRun.txt --suffix aTest --harvest
```

**c) draw the plots**
With the `dataMCcomparison.C` macro
