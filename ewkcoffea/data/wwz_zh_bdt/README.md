
## Information on TMVA to XGBoost conversion
To convert the TMVA to XGB, the [tmva-to-xgboost](https://github.com/guitargeek/tmva-to-xgboost) tool is used. On HPG, some setup is required, you can copy this into a setup shell script (or just run the commands from the command line):
```
#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

export SCRAM_ARCH=el8_amd64_gcc10
export CMSSW_VERSION=CMSSW_12_5_0
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /cvmfs/cms.cern.ch/$SCRAM_ARCH/cms/cmssw/$CMSSW_VERSION/src
eval `scramv1 runtime -sh`
cd - > /dev/null

echo 'Setup following ROOT'
which root

#eof
```
To compile, clone the repo and cd into it then:
```
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/cvmfs/cms.cern.ch/el8_amd64_gcc10/external/boost/1.80.0-2266b8f9c5c3ec1784bdd015fb86b58f/lib/
g++ -o tmva2xgboost -g -pthread -std=c++17 -m64 -I/cvmfs/cms.cern.ch/el8_amd64_gcc10/external/boost/1.80.0-2266b8f9c5c3ec1784bdd015fb86b58f/include/ -L/cvmfs/cms.cern.ch/el8_amd64_gcc10/external/boost/1.80.0-2266b8f9c5c3ec1784bdd015fb86b58f/lib/ -lboost_program_options -L/cvmfs/cms.cern.ch/el8_amd64_gcc10/external/boost/1.80.0-2266b8f9c5c3ec1784bdd015fb86b58f/lib/ tmva2xgboost.cpp
```
Next you can run the `tmva2xgboost` as described in the tool's readme, e.g.:
```
./tmva2xgboost --input sf_wwz_wgt.xml --n_features 26 > test.json
```

## For the SR BDTs 

This is where the current (v5) BDT xml files were downloaded from: 
```
#OF: 
wget http://uaf-10.t2.ucsd.edu/~kdownham/WWZ/BDT_Information/v5_020824/weights/Opposite_Flavor/WWZ/020824_newBkgd_BDTG_LR0p1_D2_NS1.weights.xml -O of_wwz_wgt.xml
wget http://uaf-10.t2.ucsd.edu/~kdownham/WWZ/BDT_Information/v5_020824/weights/Opposite_Flavor/ZH/020824_newBkgd_BDTG_T400_LR0p5.weights.xml    -O of_zh_wgt.xml

#SF: 
wget http://uaf-10.t2.ucsd.edu/~kdownham/WWZ/BDT_Information/v5_020824/weights/Same_Flavor/WWZ/020824_newBkgd_BDTG_LR0p1.weights.xml  -O sf_wwz_wgt.xml
wget http://uaf-10.t2.ucsd.edu/~kdownham/WWZ/BDT_Information/v5_020824/weights/Same_Flavor/ZH/020824_newBkgd_BDTG_LR0p1.weights.xml   -O sf_zh_wgt.xml
```
Then using the TMVA to XGBoost setup described above, the following commands were run (note if number of features changes in future versions of the model, will need to specify different n_features):
```
./tmva2xgboost --input sf_wwz_wgt.xml --n_features 27 > sf_WWZ.json
./tmva2xgboost --input sf_zh_wgt.xml  --n_features 27 > sf_ZH.json
./tmva2xgboost --input of_wwz_wgt.xml --n_features 26 > of_WWZ.json
./tmva2xgboost --input of_zh_wgt.xml  --n_features 27 > of_ZH.json
```
