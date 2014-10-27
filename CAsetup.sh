#!/bin/sh
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
cd Validation/RecoTrack/test
cp $DIR/Validation/RecoTrack/test/dummyTestLight_cfg.py Validation/RecoTrack/test/
cp $DIR/Validation/RecoTrack/test/MultiTrackValidatorForCA_cfg.py Validation/RecoTrack/test/
echo 'make some changes in SeedingHitSet.h'
cp $DIR/RecoTracker/TkSeedingLayers/interface/SeedingHitSet.h RecoTracker/TkSeedingLayers/interface/
echo 'make some changes in HitExtractorSTRP'
cp $DIR/RecoTracker/TkSeedingLayers/src/HitExtractorSTRP.h /RecoTracker/TkSeedingLayers/src/
cp $DIR/RecoTracker/TkSeedingLayers/src/HitExtractorSTRP.cc /RecoTracker/TkSeedingLayers/src/

echo 'Copying CA code'
cp -r $DIR/RecoTracker/CAtracker/ RecoTracker/CAtracker/

echo 'building the code'
scram b
echo 'done! To make a simple test go to  Validation/RecoTrack/test and run: cmsRun dummyTestLight_cfg.py'
