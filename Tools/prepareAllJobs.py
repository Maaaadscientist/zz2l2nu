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

def copy_catalog_files_on_local(theCatalog, jobID, jobSpliting):
    scriptLines = ''
    try:
        theCatalogFile = open(theCatalog,'r')
    except IOError:
        sys.stderr.write("cannot open catalog file")

    iteFileInJob=0
    minIteJob = jobID*jobSpliting
    maxIteJob = (jobID+1)*jobSpliting

    theCatalogLines = theCatalogFile.readlines()
    for aLine in theCatalogLines:
        if not ".root" in aLine: continue
        if (aLine.startswith("#")): continue  
        iteFileInJob=iteFileInJob+1
        if (iteFileInJob<=minIteJob) or (iteFileInJob>maxIteJob): continue
        fileName = re.split(" ",aLine)[0]
        scriptLines += 'dccp '+fileName+' .\n'

    scriptLines += 'ls *.root > theLocalCata.txt\n'
    return scriptLines

def extract_list_of_systs(syst):
    dictOfSysts = {}
    if not syst or syst=="no":
      dictOfSysts[None] = [""]
    elif ("_up" in syst) or ("_down" in syst):
      dictOfSysts[syst] = [""]
    elif syst=="all":
      dictOfSysts=parse_syst_file()
      dictOfSysts[None] = [""]
    else:
      dictOfSysts[syst+"_up"] = [""]
      dictOfSysts[syst+"_down"] = [""]
    return dictOfSysts

def prepare_job_script(theCatalog, name,jobID,isMC,jobSpliting,currentSyst):
    global base_path
    global thisSubmissionDirectory
    global outputDirectory
    global jobsDirectory

    scriptFile = open(jobsDirectory+'/scripts/runOnBatch_'+outputPrefixName+name+'_'+str(jobID)+'.sh','w')
    scriptLines = ''
    scriptLines += 'source $VO_CMS_SW_DIR/cmsset_default.sh\n'
    scriptLines += 'export SCRAM_ARCH=slc6_amd64_gcc530\n'
#    scriptLines += 'export BUILD_ARCH=slc6_amd64_gcc530\n'
#    scriptLines += 'export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch\n'
#    scriptLines += 'export XRD_NETWORKSTACK=IPv4\n'
    scriptLines += ('export INITDIR='+base_path+'\n')
    scriptLines += ('cd $INITDIR\n')
    scriptLines += 'ulimit -s 20480\n'
    scriptLines += 'eval `scramv1 runtime -sh`\n'
    scriptLines += 'cd -\n'
#    scriptLines += 'ulimit -c 0;\n'
    scriptLines += 'if [ -d $TMPDIR ] ; then cd $TMPDIR ; fi;\n'
    scriptLines += 'cp '+thisSubmissionDirectory+'/runHZZanalysis .;\n'
    scriptLines += 'cp -r '+base_path+'/data .;\n'
    scriptLines += 'hostname ;\n'
    #iteFileInJob=0
#    for aFile in listFiles:
    scriptLines += ("date;\n")
#        scriptLines += ("dccp "+aFile+" inputFile_"+str(jobID)+"_"+str(iteFileInJob)+".root;\n")
    keepAllControlPlotsOption = ""
    if args.syst=="all":
      keepAllControlPlotsOption = " keepAllControlPlots=false"
    if args.localCopy:
        scriptLines += copy_catalog_files_on_local(theCatalog, jobID, jobSpliting)
        if currentSyst:
          scriptLines += ("./runHZZanalysis catalogInputFile=theLocalCata.txt histosOutputFile="+outputPrefixName+name+"_"+str(jobID)+".root skip-files=0 max-files="+str(jobSpliting)+" isMC="+str(isMC)+" maxEvents=-1 isPhotonDatadriven="+str(isPhotonDatadriven)+" doInstrMETAnalysis="+str(doInstrMETAnalysis)+" doTnPTree="+str(doTnPTree)+" syst="+currentSyst+keepAllControlPlotsOption+";\n")
        else:
          scriptLines += ("./runHZZanalysis catalogInputFile=theLocalCata.txt histosOutputFile="+outputPrefixName+name+"_"+str(jobID)+".root skip-files=0 max-files="+str(jobSpliting)+" isMC="+str(isMC)+" maxEvents=-1 isPhotonDatadriven="+str(isPhotonDatadriven)+" doInstrMETAnalysis="+str(doInstrMETAnalysis)+" doTnPTree="+str(doTnPTree)+keepAllControlPlotsOption+";\n")
    else:
        if currentSyst:
          scriptLines += ("./runHZZanalysis catalogInputFile="+theCatalog+" histosOutputFile="+outputPrefixName+name+"_"+str(jobID)+".root skip-files="+str(jobID*jobSpliting)+" max-files="+str(jobSpliting)+" isMC="+str(isMC)+" maxEvents=-1 isPhotonDatadriven="+str(isPhotonDatadriven)+" doInstrMETAnalysis="+str(doInstrMETAnalysis)+" doTnPTree="+str(doTnPTree)+" syst="+currentSyst+keepAllControlPlotsOption+";\n")
        else:
          scriptLines += ("./runHZZanalysis catalogInputFile="+theCatalog+" histosOutputFile="+outputPrefixName+name+"_"+str(jobID)+".root skip-files="+str(jobID*jobSpliting)+" max-files="+str(jobSpliting)+" isMC="+str(isMC)+" maxEvents=-1 isPhotonDatadriven="+str(isPhotonDatadriven)+" doInstrMETAnalysis="+str(doInstrMETAnalysis)+" doTnPTree="+str(doTnPTree)+keepAllControlPlotsOption+";\n")
#        scriptLines += ("rm inputFile_"+str(jobID)+"_"+str(iteFileInJob)+".root;\n\n")
#        iteFileInJob = iteFileInJob+1
#    scriptLines += ('$ROOTSYS/bin/hadd output_'+name+"_"+str(jobID)+".root theOutput_"+name+"_"+str(jobID)+"_*.root;\n\n")
    scriptLines += ("cp "+outputPrefixName+name+"_"+str(jobID)+".root "+outputDirectory+"\n")
    scriptFile.write(scriptLines)
    scriptFile.close()

    #jobsFiles = open("sendJobs_"+re.split("_",outputDirectory)[1]+".cmd","a")
    jobsFiles = open(thisSubmissionDirectory+"/sendJobs_"+args.suffix+".cmd","a")
    jobsFiles.write("qsub "+str(doExpress)+" -j oe -o "+jobsDirectory+'/logs/ '+jobsDirectory+'/scripts/runOnBatch_'+outputPrefixName+name+'_'+str(jobID)+'.sh\n')
    jobsFiles.close()
    #print scriptLines

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
    jobSpliting=25
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
            curentSize = curentSize+int(lineField[1])
            if len(listFileInAJob)>=jobSpliting: #curentSize>5000000000:
                #print("jobID="+str(jobID))
                prepare_job_script(catalogDirectory+'/'+catalogName, shortName+systString, jobID, isMC, jobSpliting, currentSyst)
                listFileInAJob=[]
                curentSize=0
                jobID+=1
    if len(listFileInAJob)>0 :
        #there are remaining files to run
        #print("jobIDr="+str(jobID))
        prepare_job_script(catalogDirectory+'/'+catalogName, shortName+systString, jobID, isMC, jobSpliting, currentSyst)


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
    global doExpress
    #create the directories if needed
    base_path=os.path.expandvars('$CMSSW_BASE/src/shears/HZZ2l2nu')
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
    shutil.copy2(base_path+'/runHZZanalysis', thisSubmissionDirectory)

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
