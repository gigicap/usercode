#!/bin/bash
echo 'this script works for CMSSW_7_4_0_pre1' 

cmsenv
#these lines will be removed in the final script
echo 'Initialising git' 
eval `ssh-agent`
ssh-add
git cms-init

echo 'adding pkg Validation/RecoTrack/test' 
git cms-addpkg Validation/RecoTrack/test
echo 'adding pkg RecoTracker/TkSeedingLayers/'
git cms-addpkg RecoTracker/TkSeedingLayers/
cd RecoTracker/TkSeedingLayers/interface
rm SeedingHitSet.h
cd ../src/
rm HitExtractorSTRP.h
rm HitExtractorSTRP.cc
cd ../../..

#DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo 'add some CA related validation scripts'
cp ../../CAcode/src/Validation/RecoTrack/test/dummyTestLight_cfg.py Validation/RecoTrack/test/dummyTestLight_cfg.py
cp ../../CAcode/src/Validation/RecoTrack/test/MultiTrackValidatorForCA_cfg.py Validation/RecoTrack/test/MultiTrackValidatorForCA_cfg.py
echo 'make some changes in SeedingHitSet.h'
cp ../../CAcode/src/RecoTracker/TkSeedingLayers/interface/SeedingHitSet.h RecoTracker/TkSeedingLayers/interface/SeedingHitSet.h
echo 'make some changes in HitExtractorSTRP'
cp ../../CAcode/src/RecoTracker/TkSeedingLayers/src/HitExtractorSTRP.h RecoTracker/TkSeedingLayers/src/HitExtractorSTRP.h
cp ../../CAcode/src/RecoTracker/TkSeedingLayers/src/HitExtractorSTRP.cc RecoTracker/TkSeedingLayers/src/HitExtractorSTRP.cc

echo 'Copying CA code'
cp -r ../../CAcode/src/RecoTracker/CAtracker/ RecoTracker/CAtracker/

echo 'Checking deps'
git cms-checkdeps -a

echo 'building the code (might take a while)'
scram b
echo 'done! To make a simple test go to  Validation/RecoTrack/test and run: cmsRun dummyTestLight_cfg.py'
