# Consistent LCG environment (http://lcginfo.cern.ch)
. /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_96python3 x86_64-slc6-gcc8-opt

export HZZ2L2NU_BASE=$(pwd)
export PYTHONPATH=$HZZ2L2NU_BASE/python:$PYTHONPATH
export PATH=$HZZ2L2NU_BASE/bin:$PATH

