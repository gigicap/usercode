Cellular Automaton code +
changes needed to make it run in CMSSW_7_4_0_pre1 +
sh script to copy everything in the right place.

========

How to run the code:

mkdir CAtracking
cd CAtracking
git clone -b CAforRelease_7_4_0 https://github.com/gigicap/usercode.git CAcode/
cmsrel CMSSW_7_4_0_pre1
cd CMSSW_7_4_0_pre8/src
cmsenv
source ../../CAcode/CAsetup.sh

the CAcode folder can then be removed.

========
