#!/bin/bash
echo 'this script works for CMSSW_7_2_0_pre8' 

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

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo 'add some CA related validation scripts'
cp $DIR/src/Validation/RecoTrack/test/dummyTestLight_cfg.py Validation/RecoTrack/test/dummyTestLight_cfg.py
cp $DIR/src/Validation/RecoTrack/test/MultiTrackValidatorForCA_cfg.py Validation/RecoTrack/test/MultiTrackValidatorForCA_cfg.py
echo 'make some changes in SeedingHitSet.h'
cp $DIR/src/RecoTracker/TkSeedingLayers/interface/SeedingHitSet.h RecoTracker/TkSeedingLayers/interface/SeedingHitSet.h
echo 'make some changes in HitExtractorSTRP'
cp $DIR/src/RecoTracker/TkSeedingLayers/src/HitExtractorSTRP.h /RecoTracker/TkSeedingLayers/src/HitExtractorSTRP.h
cp $DIR/src/RecoTracker/TkSeedingLayers/src/HitExtractorSTRP.cc /RecoTracker/TkSeedingLayers/src/HitExtractorSTRP.cc

echo 'Copying CA code'
cp -r $DIR/src/RecoTracker/CAtracker/ RecoTracker/CAtracker/

echo 'building the code'
scram b
echo 'done! To make a simple test go to  Validation/RecoTrack/test and run: cmsRun dummyTestLight_cfg.py'
