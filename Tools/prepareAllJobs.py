#!/usr/bin/env python



import sys
import re
import os
import argparse
import shutil

def parse_command_line():
    """Parse the command line parser"""
    parser = argparse.ArgumentParser(description='Launch baobab nutple production.')

    parser.add_argument('--listDataset', action='store', default=None,
                        help='Specifies the file containing the list of task to submit. The file must contain one task specification per line. The specification consist of three space-separated value: catalog of input files, pruner selection, pruner subselection. The catalog can be an eos path (/store/...).')
    parser.add_argument('--suffix', action='store', default=None,
                        help='suffix that will be added to the output directory')
    parser.add_argument('--harvest', action='store_true', default=None,
                        help='harvest the root files from the last submission')
    parser.add_argument('--isPhotonDatadriven', action='store_true', default=None,
                        help='Launch HZZ with Instr. MET DY estimated from photon')
    parser.add_argument('--doInstrMETAnalysis', action='store_true', default=None,
                        help='Launch InstrMETAnalysis')
    parser.add_argument('--doTnPTree', action='store_true', default=None,
                        help='Launch TnP Tree production')
    parser.add_argument('--doNRBAnalysis', action='store_true', default=None,
                        help='Launch NRB Analysis')
    parser.add_argument('--express', action='store_true', default=None,
        help='Launch on the express queue (better if: very fast jobs "<10min" in a small amount "<16"). This queue has a wallTime of 32min and can only take 8 jobs per user')
    parser.add_argument('--localCopy', action='store_true', default=None,
        help='Copy the ROOT files locally on the node to improve speed.')
    parser.add_argument('--syst', action='store', default=None,
        help='Specify the systematic on which you need to run. If you dont specify _up or _down, both will be run at the same time. Use "all" to run on all the systematics defined in the systList.txt file.')

    return parser.parse_args()

def parse_datasets_file():
    global catalogDirectory
    try:
    	datasetFile = open(args.listDataset,'r')
    except KeyError:
        sys.stderr.write("please specify a list of datasets")

    datasets = datasetFile.readlines()
    listCatalogs=[]
    for aLine in datasets:
        if (aLine.startswith("#")): #allow comments in the dataset list
            print("\033[1;32m Ignoring line \033[1;34m"+aLine[:-1]+"\033[0;m")
        elif ("catalogPath=" in aLine):
            catalogDirectory = (re.split("=",aLine)[1])[:-1]
            print("\033[1;32m the catalogs are in the directory \033[1;34m"+catalogDirectory+"\033[0;m")
        else:
            listCatalogs.append(aLine[:-1])
    return listCatalogs

def parse_syst_file():
    global base_path
    try:
      systFile = open(base_path+"/systList.txt")
    except KeyError:
      sys.stderr.write("cannot open syst file")
    systLines = systFile.readlines()
    theSystDict={}
    for aLine in systLines:
      if aLine.startswith("//"): continue
      if aLine.startswith("#"): continue
      systKey=(aLine.split(" ")[0])
      thisLineAsList = aLine.split()
      theSystDict[systKey] = []
      for i in thisLineAsList:
        if i==thisLineAsList[0]: continue
        if i.startswith("//"): break
        if i.startswith("#"): break
        theSystDict[systKey].append(i)
    return theSystDict

def find_syst_in_file(currentSyst):
    global base_path
    try:
      systFile = open(base_path+"/systList.txt")
    except KeyError:
      sys.stderr.write("cannot open syst file")
    systLines = systFile.readlines()
    thisSystList=[]
    for aLine in systLines:
      if aLine.startswith("//"): continue
      if aLine.startswith("#"): continue
      if not currentSyst in aLine: continue
      thisLineAsList = aLine.split()
      for i in thisLineAsList:
        if i==thisLineAsList[0]: continue
        if i.startswith("//"): break
        if i.startswith("#"): break
        thisSystList.append(i)
    return thisSystList
    

def copy_catalog_files_on_local(catalog_path, job_id, job_splitting):
    """Construct shell commands to copy a subset of files from a catalog

    Read paths of files in the given catalog and select those that
    correspond to the given job.  Construct shell commands to copy the
    selected files to the current directory.

    Arguments:
        catalog_path:  Path to a catalog file.
        job_id:  Zero-based index of the requested job.
        job_splitting:  Number of files to be processed per job.

    Return value:
        List of strings with the copy commands as well as a command to
        create a catalog of copied files.
    """

    try:
        with open(catalog_path, 'r') as f:
            catalog_lines = f.readlines()
    except IOError:
        sys.stderr.write('Cannot open catalog file "{}".'.format(catalog_path))

    catalog_files = []

    for line in catalog_lines:
        if line.startswith('#') or '.root' not in line:
            continue

        if 'Bonzai' in line:
            file_name = line.split()[0]
        else:
            file_name = line.strip()

        catalog_files.append(file_name)

    script_commands = []

    for i in range(
        job_id * job_splitting,
        min((job_id + 1) * job_splitting, len(catalog_files))
    ):
        script_commands.append('dccp {} .'.format(catalog_files[i]))

    script_commands.append('ls *.root > theLocalData.txt')
    return script_commands


def extract_list_of_systs(syst):
    dictOfSysts = {}
    if not syst or syst=="no":
      dictOfSysts[None] = [""]
    elif ("_up" in syst) or ("_down" in syst):
      dictOfSysts[syst] = find_syst_in_file(syst)
    elif syst=="all":
      dictOfSysts=parse_syst_file()
      dictOfSysts[None] = [""]
    else:
      dictOfSysts[syst+"_up"] = find_syst_in_file(syst+"_up")
      dictOfSysts[syst+"_down"] = find_syst_in_file(syst+"_down")
      dictOfSysts[None] = [""] # Run also on the nominal shape
    return dictOfSysts


def prepare_job_script(
    catalog_path, name, job_id, is_mc, job_splitting, current_syst
):
    """Create a script to be run on a batch system

    The script performs set-up and executes runHZZanalysis.  It is saved
    in the jobs directory.

    Arguments:
        catalog_path:  Path to a catalog file.
        name:    Name for the task.
        job_id:  Zero-based index of requested job.
        is_mc:   Indicates whether simulation or real data are being
            processed.
        job_splitting:  Number of files to process per job.
        current_syst:   Label of requested systematic variation.

    Return value:
        None.
    """

    global base_path
    global thisSubmissionDirectory
    global outputDirectory
    global jobsDirectory

    script_commands = [
      'export INITDIR={}\n'.format(base_path),
      'cd $INITDIR',
      '. ./env.sh',
      'cd -',
      # 'ulimit -c 0',
      'if [ -d $TMPDIR ] ; then cd $TMPDIR ; fi',
      'cp {}/runHZZanalysis .'.format(thisSubmissionDirectory),
      'cp -r {}/data .'.format(base_path),
      'hostname',
      'date'
    ]

    if args.localCopy:
        script_commands += copy_catalog_files_on_local(
            catalog_path, job_id, job_splitting
        )

    # Construct options for runHZZanalysis program
    options = [
        '--catalog={}'.format(
            'theLocalData.txt' if args.localCopy else catalog_path
        ),
        '--output={}{}_{}.root'.format(outputPrefixName, name, job_id),
        '--skip-files={}'.format(
            0 if args.localCopy else job_id * job_splitting
        ),
        '--max-files={}'.format(job_splitting), '--max-events=-1',
        '--is-mc={}'.format(is_mc)
    ]

    if doInstrMETAnalysis:
        options.append('--analysis=InstrMET')
    elif doTnPTree:
        options.append('--analysis=TnP')
    elif doNRBAnalysis:
        options.append('--analysis=NRB')
    else:
        options.append('--analysis=Main')

    if isPhotonDatadriven:
        options.append('--dd-photon')

    if current_syst:
        options.append('--syst={}'.format(current_syst))

    if args.syst != 'all':
        options.append('--all-control-plots')

    # Use debug-level verbosity when running on the batch system since
    # one is expected to consult logs only in case of a problem
    options.extend(['-v', '2'])

    run_application_command = ' '.join(['./runHZZanalysis'] + options)
    script_commands.append('echo ' + run_application_command)
    script_commands.append(run_application_command)

    # script_commands.append(
    #     '$ROOTSYS/bin/hadd output_{name}_{jobid}.root '
    #     'theOutput_{name}_{jobid}_*.root'.format(name=name, jobid=job_id)
    # )
    script_commands.append(
        'cp {}{}_{}.root {}'.format(outputPrefixName, name, job_id, outputDirectory)
    )


    script_path = '{}/scripts/runOnBatch_{}{}_{}.sh'.format(
        jobsDirectory, outputPrefixName, name, job_id
    )

    with open(script_path, 'w') as f:
        for command in script_commands:
            f.write(command)
            f.write('\n')

    with open(
        '{}/sendJobs_{}.cmd'.format(thisSubmissionDirectory, args.suffix), 'a'
    ) as jobs_file:
        jobs_file.write(
            'qsub {} -j oe -o {}/logs/ {}\n'.format(
                doExpress, jobsDirectory, script_path
            )
        )


def make_the_name_short(theLongName):
    shortName=''
    shortNameIntermediate = ((theLongName.split("-", 1)[1]).split("Pruner")[0]).rsplit("-", 1)[0]
    shortName = (shortNameIntermediate.split("_TuneCUETP8M1")[0]).split("_13TeV")[0]
    return shortName

def create_script_fromCatalog(catalogName,currentSyst):
    isMC = 0
    shortName=make_the_name_short(catalogName)

    print("Preparing the scripts for \033[1;33m"+catalogName+"\033[0;m with short name=\033[1;33m"+shortName+"\033[0;m")

    catalogFile = open(catalogDirectory+'/'+catalogName,'r')
    catalogLines = catalogFile.readlines()

    curentSize = 0
    listFileInAJob=[]
    jobID=0
    jobSplitting=25
    if 'Baobab' in catalogName:
        jobSplitting=10
    if (currentSyst and ("pdf" in currentSyst or "QCDscale" in currentSyst)): jobSplitting=9999
    jobID=0
    listFileInAJob=[]
    curentSize=0
    if not currentSyst:
      systString = ""
    else:
      systString = '_'+currentSyst
    for aLine in catalogLines:
        if ("data type" in aLine):
            if ("mc" in aLine):
                #print("this sample is a MC sample")
                isMC = 1
        if ".root" in aLine:
            lineField=re.split(" ",aLine)
            listFileInAJob.append(lineField[0])
            if "bonzai" in aLine:
                curentSize = curentSize+int(lineField[1])
            else:
                curentSize = curentSize+200000000
            if len(listFileInAJob)>=jobSplitting: #curentSize>5000000000:
                #print("jobID="+str(jobID))
                prepare_job_script(catalogDirectory+'/'+catalogName, shortName+systString, jobID, isMC, jobSplitting, currentSyst)
                listFileInAJob=[]
                curentSize=0
                jobID+=1
    if len(listFileInAJob)>0 :
        #there are remaining files to run
        #print("jobIDr="+str(jobID))
        prepare_job_script(catalogDirectory+'/'+catalogName, shortName+systString, jobID, isMC, jobSplitting, currentSyst)


def runHarvesting():
    global thisSubmissionDirectory
    global outputDirectory

    try:
        datasetFile = open(args.listDataset,'r')
    except KeyError:
        sys.stderr.write("please specify a list of datasets")
    if not os.path.isdir(thisSubmissionDirectory+"/MERGED"):
      print("\033[1;34m will create the directory "+thisSubmissionDirectory+"/MERGED"+"\033[0;m")
      os.mkdir(thisSubmissionDirectory+"/MERGED")
    dataSamplesList = ""
    listForFinalPlots = {}
    listForFinalPlots_data = ""
    dictOfSysts = extract_list_of_systs(args.syst)
    for currentSyst in dictOfSysts:
      datasetFile.seek(0)
      dataSamplesList = ""
      dataForThisSyst = None
      if not currentSyst:
        systString = ""
      else:
        systString = '_'+currentSyst
      for aLine in datasetFile:
        harvestForThisSyst = None
        for key in dictOfSysts[currentSyst]:
          if key in aLine: harvestForThisSyst = True
        if not harvestForThisSyst: continue
        if (aLine.startswith("#")): continue
        if (aLine.startswith("catalogPath")): continue
        if not "Bonzais" in aLine: continue
        theShortName=make_the_name_short(aLine[:-1])
        print("\033[1;32m merging "+theShortName+systString+"\033[0;m")
        os.system("$ROOTSYS/bin/hadd -f "+thisSubmissionDirectory+"/MERGED/"+outputPrefixName+theShortName+systString+".root "+outputDirectory+"/"+outputPrefixName+theShortName+systString+"_[0-9]*.root")
        if "Data" in aLine:
          dataForThisSyst = True
          dataSamplesList = dataSamplesList+" "+thisSubmissionDirectory+"/MERGED/"+outputPrefixName+theShortName+systString+".root"
        else:
          if theShortName in listForFinalPlots:
            listForFinalPlots[theShortName] = listForFinalPlots[theShortName]+" "+thisSubmissionDirectory+"/MERGED/"+outputPrefixName+theShortName+systString+".root"
          else:
            listForFinalPlots[theShortName] = thisSubmissionDirectory+"/MERGED/"+outputPrefixName+theShortName+systString+".root"
      if dataForThisSyst: listForFinalPlots_data = listForFinalPlots_data + " "+thisSubmissionDirectory+"/MERGED/"+outputPrefixName+"Data"+systString+".root"
      if currentSyst: print("\033[1;32m merging all Data (Single* and Double*) together for "+currentSyst+"\033[0;m")
      else: print("\033[1;32m merging all Data (Single* and Double*) together for nominal shapes\033[0;m")
      if dataForThisSyst: os.system("$ROOTSYS/bin/hadd -f "+thisSubmissionDirectory+"/MERGED/"+outputPrefixName+"Data"+systString+".root "+dataSamplesList)
    if args.syst and not args.syst=="no":
      print("\033[1;32m producing final ROOT files with all shapes\033[0;m")
      for key in listForFinalPlots:
        os.system("$ROOTSYS/bin/hadd -f "+thisSubmissionDirectory+"/MERGED/"+outputPrefixName+key+"_final.root "+listForFinalPlots[key])
      os.system("$ROOTSYS/bin/hadd -f "+thisSubmissionDirectory+"/MERGED/"+outputPrefixName+"Data_final.root "+listForFinalPlots_data)


def main():
    global args
    global catalogDirectory
    global base_path
    global thisSubmissionDirectory
    global outputDirectory
    global jobsDirectory
    global doInstrMETAnalysis
    global isPhotonDatadriven
    global outputPrefixName
    global doTnPTree
    global doNRBAnalysis
    global doExpress
    #create the directories if needed
    base_path=os.path.expandvars('$HZZ2L2NU_BASE')
    if not os.path.isdir(base_path+"/OUTPUTS"):
        print("\033[1;31m OUTPUTS directory does not exist: will create it \033[0;m")
        os.mkdir(base_path+"/OUTPUTS")

    args = parse_command_line()
    
    if type(args.suffix) != type("txt"):
        thisSubmissionDirectory=base_path+"/OUTPUTS/Test"
    else:
        thisSubmissionDirectory=base_path+"/OUTPUTS/"+args.suffix
    
    outputDirectory=thisSubmissionDirectory+"/OUTPUTS"
    jobsDirectory=thisSubmissionDirectory+"/JOBS"
    
    if not os.path.isdir(thisSubmissionDirectory):
        print("\033[1;31m will create the directory "+thisSubmissionDirectory+"\033[0;m")
        os.mkdir(thisSubmissionDirectory)
    if not os.path.isdir(outputDirectory):
        os.mkdir(outputDirectory)
    if not os.path.isdir(jobsDirectory):
        os.mkdir(jobsDirectory)
        os.mkdir(jobsDirectory+"/scripts")
        os.mkdir(jobsDirectory+"/logs")

    #options
    outputPrefixName="outputHZZ_"
    if args.doInstrMETAnalysis:
        print "Preparing InstrMET analysis...\n"
        doInstrMETAnalysis = 1
        outputPrefixName="outputInstrMET_"
    else:
        doInstrMETAnalysis = 0

    if args.isPhotonDatadriven:
        print "Datadriven estimation of the Instr. MET option found...\n"
        isPhotonDatadriven = 1
        outputPrefixName="outputPhotonDatadriven_"
    else:
        isPhotonDatadriven = 0


    if args.doTnPTree:
        print "Praparing Tag and Probe Tree...\n"
        doTnPTree = 1
        outputPrefixName="outputTnP_"
    else:
        doTnPTree = 0
    
    if args.doNRBAnalysis:
        print "Praparing Non-resonant Bkg. Analysis...\n"
        doNRBAnalysis = 1
        outputPrefixName="outputNRB_"
    else:
        doNRBAnalysis = 0
    if args.harvest:
        print "will harvest"
        runHarvesting()
        return

    if args.express:
        print "Will be launched on the express queue (NB: only do this for small and fast jobs)\n"
        doExpress = " -q express -l walltime=00:30:00 "
    else:
        print "WallTime is set to 20h. If you need more, please update the script. If you need to send only a small number of very short jobs, please consider using the express queue (--express)\n"
        doExpress = " -l walltime=20:00:00 "



    listCatalogs=parse_datasets_file()

    #copy catalog list and executable to the OUTPUTS directory so we can run in parallel and always have a backup of what we ran
    #shutil.copy2(args.listDataset, thisSubmissionDirectory+'/'+os.path.basename(args.listDataset)) #This is now done in the launchAnalysis script
    shutil.copy2(base_path+'/bin/runHZZanalysis', thisSubmissionDirectory)

    #check if the file for big submission does exist and then remove it
    #Hugo: the way the Instr. MET is done, I'm updating the big submission script so please don't remove it while preparing jobs.
    #if os.path.exists("sendJobs_"+re.split("_",outputDirectory)[1]+".cmd"):
    #    print("\033[1;31m sendJobs_"+re.split("_",outputDirectory)[1]+".cmd already exist-> removing it ! \033[0;37m")
    #    os.remove("sendJobs_"+re.split("_",outputDirectory)[1]+".cmd")

    dictOfSysts = extract_list_of_systs(args.syst)

    for thisSyst in dictOfSysts:
      for aCatalog in listCatalogs:
        runOnThisCatalog = None
        for key in dictOfSysts[thisSyst]:
          if key in aCatalog: runOnThisCatalog = True
        if runOnThisCatalog: create_script_fromCatalog(aCatalog,thisSyst)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt, e:
        print "\nBye!"
        pass
