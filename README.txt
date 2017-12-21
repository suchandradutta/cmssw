# Print DetIds 
cd SLHCUpgradeSimulations/Geometry/test 
cmsRun printGeometry_pixelTelescope_cfg.py


# Create SimHitMap
# Needs digi geometry !! (not there yet)
cmsDriver.py MinBias_TuneZ2star_14TeV_pythia6_cff --geometry Phase2TestBeam --conditions auto:phase2_realistic -n 200 --era Phase2 --eventcontent FEVTDEBUG --relval 10000,100 -s GEN,SIM --datatier GEN-SIM --beamspot HLLHC --no_exec
gedit MinBias_TuneZ2star_14TeV_pythia6_cff_GEN_SIM.py &   #edit path to reco geom file
# Can only run that if no change in geometry cff file
cmsRun MinBias_TuneZ2star_14TeV_pythia6_cff_GEN_SIM.py


cmsRun Demo/SimHitAnalyzer/test/TTbar_cfi_GEN_SIM.py
root hsimhit_minbias_200.root
