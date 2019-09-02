# ShowerShapeTest


```
cmsrel CMSSW_10_2_13 #or any modern CMSSW release should work in most
cd CMSSW_10_2_13/src
git clone git@github.com:Sam-Harper/ShowerShapeTest.git 
cd ShowerShapeTest
scram b -j 16
cmsRun ShowerShapeTest/test/runLocalEnergyMapProducer.py inputFiles=input.root outputFile=output.root
```
