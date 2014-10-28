Cellular Automaton code +
changes needed to make it run in CMSSW_7_2_0_pre8 +
sh script to copy everything in the right place.

========

How to run the code:

mkdir CAtracking
cd CAtracking
git clone -b master_simplified https://github.com/gigicap/usercode.git CAcode/
cmsrel CMSSW_7_2_0_pre8
cd CMSSW_7_2_0_pre8/src
cmsenv
source ../../CAcode/CAsetup.sh

the CAcode folder can then be removed.

========
