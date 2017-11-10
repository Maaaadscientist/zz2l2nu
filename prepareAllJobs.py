#!/usr/bin/env python



import sys
import re
import os
import argparse


def parse_command_line():
    """Parse the command line parser"""
    parser = argparse.ArgumentParser(description='Launch baobab nutple production.')

    parser.add_argument('--listDataset', action='store', default=None,
                        help='Specifies the file containing the list of task to submit. The file must contain one task specification per line. The specificaiton consist of three space-separated value: catalog of input files, pruner selection, pruner subselection. The catalog can be an eos path (/store/...).')
    parser.add_argument('--suffix', action='store', default=None,
                        help='suffix that will be added to the output directory')
    parser.add_argument('--harvest', action='store_true', default=None,
                        help='harvest the root files from the last submission')

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
        if ("catalogPath=" in aLine):
            catalogDirectory = (re.split("=",aLine)[1])[:-1]
            print("\033[1;32m the catalogs are in the directory \033[1;34m"+catalogDirectory+"\033[1;37m")
        else:
            listCatalogs.append(aLine[:-1])
    return listCatalogs


def prepare_job_script(theCatalog, name,jobID,isMC,jobSpliting):
    global outputDirectory

    scriptFile = open('JOBS/runOnBatch_'+name+'_'+str(jobID)+'.sh','w')
    scriptLines = ''
    scriptLines += 'source $VO_CMS_SW_DIR/cmsset_default.sh\n'
    scriptLines += 'export SCRAM_ARCH=slc6_amd64_gcc530\n'
#    scriptLines += 'export BUILD_ARCH=slc6_amd64_gcc530\n'
#    scriptLines += 'export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch\n'
#    scriptLines += 'export XRD_NETWORKSTACK=IPv4\n'
    scriptLines += ('export INITDIR='+os.getcwd()+'\n')
    scriptLines += ('cd $INITDIR\n')
    scriptLines += 'eval `scramv1 runtime -sh`\n'
    scriptLines += 'cd -\n'
#    scriptLines += 'ulimit -c 0;\n'
    scriptLines += 'if [ -d $TMPDIR ] ; then cd $TMPDIR ; fi;\n'
    scriptLines += 'cp '+os.getcwd()+'/runHZZanalysis .;\n'
    scriptLines += 'hostname ;\n'
    iteFileInJob=0
#    for aFile in listFiles:
    scriptLines += ("date;\n")
#        scriptLines += ("dccp "+aFile+" inputFile_"+str(jobID)+"_"+str(iteFileInJob)+".root;\n")
    scriptLines += ("./runHZZanalysis catalogInputFile="+theCatalog+" histosOutputFile=output_"+name+"_"+str(jobID)+".root skip-files="+str(jobID*jobSpliting)+" max-files="+str(jobSpliting)+" isMC="+str(isMC)+" maxEvents=-1;\n")
#        scriptLines += ("rm inputFile_"+str(jobID)+"_"+str(iteFileInJob)+".root;\n\n")
#        iteFileInJob = iteFileInJob+1
#    scriptLines += ('$ROOTSYS/bin/hadd output_'+name+"_"+str(jobID)+".root theOutput_"+name+"_"+str(jobID)+"_*.root;\n\n")
    scriptLines += ("cp output_"+name+"_"+str(jobID)+".root $INITDIR/"+outputDirectory+"\n")
    scriptFile.write(scriptLines)
    scriptFile.close()

    jobsFiles = open("sendJobs_"+re.split("_",outputDirectory)[1]+".cmd","a")
    jobsFiles.write("qsub -l walltime=20:00:00 -j oe "+os.getcwd()+'/JOBS/runOnBatch_'+name+'_'+str(jobID)+'.sh\n')
    jobsFiles.close()
    #print scriptLines

def make_the_name_short(theLongName):
    shortName=''
    shortNameIntermediate = re.split("-",theLongName)[1]
    if "_" in shortNameIntermediate:
        shortName = re.split("_",shortNameIntermediate)[0]
    else:
        shortName = shortNameIntermediate
    return shortName

def create_script_fromCatalog(catalogName):
    isMC = 0
    shortName=make_the_name_short(catalogName)

    print("Preparing the scripts for \033[1;33m"+catalogName+"\033[1;37m with short name=\033[1;33m"+shortName+"\033[1;37m")

    catalogFile = open(catalogDirectory+'/'+catalogName,'r')
    catalogLines = catalogFile.readlines()

    curentSize = 0
    listFileInAJob=[]
    jobID=0
    jobSpliting=25
    for aLine in catalogLines:
        if ("data type" in aLine):
            if ("mc" in aLine):
                print("this sample is a MC sample")
                isMC = 1
        if "dcap" in aLine:
            lineField=re.split(" ",aLine)
            listFileInAJob.append(lineField[0])
            curentSize = curentSize+int(lineField[1])
            if len(listFileInAJob)>jobSpliting: #curentSize>5000000000:
                print("jobID="+str(jobID))
                prepare_job_script(catalogDirectory+'/'+catalogName, shortName, jobID, isMC, jobSpliting)
                listFileInAJob=[]
                curentSize=0
                jobID+=1
    if len(listFileInAJob)>0 :
         #there are remaining files to run
        print("jobIDr="+str(jobID))
        prepare_job_script(catalogDirectory+'/'+catalogName, shortName, jobID, isMC, jobSpliting)


def runHarvesting():
    try:
        datasetFile = open(args.listDataset,'r')
    except KeyError:
        sys.stderr.write("please specify a list of datasets")
    for aLine in datasetFile:
        if not "Bonzais" in aLine:
            continue
        theShortName=make_the_name_short(aLine[:-1])
        print("\033[1;32m merging"+theShortName+"\033[1;37m")
        if not os.path.isdir("merged_"+args.suffix):
            print("\033[1;31m will create the directory "+"merged_"+args.suffix+"\033[0;37m")
            os.mkdir("merged_"+args.suffix)
        os.system("$ROOTSYS/bin/hadd -f merged_"+args.suffix+"/output_"+theShortName+".root OUTPUTS_"+args.suffix+"/output_"+theShortName+"_*.root")



def main():
    global args
    global catalogDirectory
    global outputDirectory
    #create the JOBS directory if needed
    if not os.path.isdir("JOBS"):
        print("\033[1;31m JOBS directory does not exist: will create it \033[0;37m")
        os.mkdir("JOBS")



    args = parse_command_line()

    if args.harvest:
        print "will harvest"
        runHarvesting()
        return

    listCatalogs=parse_datasets_file()
    if type(args.suffix) != type("txt"):
        outputDirectory="OUTPUTS"
    else:
        outputDirectory="OUTPUTS_"+args.suffix
    if not os.path.isdir(outputDirectory):
        print("\033[1;31m will create the directory "+outputDirectory+"\033[0;37m")
        os.mkdir(outputDirectory)

    #check if the file for big submission does exist and then remove it
    if os.path.exists("sendJobs_"+re.split("_",outputDirectory)[1]+".cmd"):
        print("\033[1;31m sendJobs_"+re.split("_",outputDirectory)[1]+".cmd already exist-> removing it ! \033[0;37m")
        os.remove("sendJobs_"+re.split("_",outputDirectory)[1]+".cmd")

    for aCatalog in listCatalogs:
        create_script_fromCatalog(aCatalog)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt, e:
        print "\nBye!"
        pass
