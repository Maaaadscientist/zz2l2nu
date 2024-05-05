export INITDIR=/user/hanwen/vbs/hzz2l2nu
cd $INITDIR
. ./env.sh
cd -
if [ -d $TMPDIR ] ; then cd $TMPDIR ; fi
hostname
date
echo runHZZanalysis --config=2016.yaml --ddf=/pnfs/iihe/cms/store/group/HZZ2l2nu/Production/2020-11-30_2016-NanoAODv7/DDF/InstrMET/GJets_HT-100To200.yaml --output=GJets_HT-100To200_1.root --skip-files=2 --max-files=2 --max-events=-1 -a PhotonTrees -v 2
runHZZanalysis --config=2016.yaml --ddf=/pnfs/iihe/cms/store/group/HZZ2l2nu/Production/2020-11-30_2016-NanoAODv7/DDF/InstrMET/GJets_HT-100To200.yaml --output=GJets_HT-100To200_1.root --skip-files=2 --max-files=2 --max-events=-1 -a PhotonTrees -v 2 || exit $?
cp GJets_HT-100To200_1.root /storage_mnt/storage/user/hanwen/vbs/hzz2l2nu/final3/output
