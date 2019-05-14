# Consistent LCG environment (http://lcginfo.cern.ch)
. /cvmfs/sft.cern.ch/lcg/views/setupViews.sh LCG_95 x86_64-slc6-gcc8-opt

# yaml-cpp is provided by LCG but doesn't seem to be included into views.
# Provide its location manually.
export YAMLCPP_ROOT_DIR="/cvmfs/sft.cern.ch/lcg/releases/yamlcpp/0.5.1-b7362/x86_64-slc6-gcc8-opt"

export HZZ2L2NU_BASE=$(pwd)
export PATH=$HZZ2L2NU_BASE/bin:$PATH

