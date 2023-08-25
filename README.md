# ewkcoffea
This analysis repository contains scripts and tools for performing analyses associated with EWK physics within the `coffea` framework, making use of the [`topcoffea`](https://github.com/TopEFT/topcoffea) package. 

## Setup instructions

First, clone the repository and `cd` into the toplevel directory. 
```
git clone https://github.com/kmohrman/ewkcoffea.git
cd ewkcoffea
```
Next, create a `conda` environment and activate it. 
```
conda env create -f environment.yml
conda activate coffea-env
```
Now we can install the `ewkcoffea` package into our new conda environment. This command should be run from the toplevel `ewkcoffea` directory, i.e. the directory which contains the `setup.py` script. 
```
pip install -e .
```
Two of the packages this analysis depends on are not conda installed (i.e. they were not included in the `environment.yml` where most of the dependencies were specified), so we can go ahead and install those into our new `conda` environment via `pip`. 
```
pip install xgboost
pip install mt2
```
The `topcoffea` package upon which this analysis also depends is not yet available on `PyPI`, so we need to clone the `topcoffea` repo and install it ourselves.
```
cd /your/favorite/directory
git clone https://github.com/TopEFT/topcoffea.git
cd topcoffea
pip install -e .  
```
Now all of the dependencies have been installed and the `ewkcoffea` repository is ready to be used. The next time you want to use it, all you have to do is to activate the environment via `conda activate coffea-env`. 

## Running the processor 

The core functionality of this analysis repository is to process NanoAOD formatted data and perform selection and produce output histograms. Assuming you are starting from the toplevel `ewkcoffea` directory and your conda environment has been activated, you can run the processor as in the following example: 
```
cd analysis/wwz/
wget -nc http://uaf-10.t2.ucsd.edu/~kmohrman/for_ci/for_wwz/WWZJetsTo4L2Nu_4F_TuneCP5_13TeV-amcatnlo-pythia8_RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2_NANOAODSIM_3LepTau_4Lep/output_1.root
python run_wwz4l.py ../../input_samples/sample_jsons/test_samples/UL17_WWZJetsTo4L2Nu_forCI.json -x futures 
``` 
This example will process a single root file (the `output_1.root` file that you just downloaded) locally using the `futures` executor. To make use of distributed resources, the `work queue` executor can be used. To use the work queue executor, just change the executor option to  `-x work_queue` and run the run script as before. Next, you will need to request some workers to execute the tasks on the distributed resources. The command depends on the batch system (and also depends on how much memory etc. you need for your particular tasks), but an example command for  `slurm` and for `condor` is given below. Note that is important to submit these commands from the same environment that you ran the manager process (so activate your  `coffea-env` in your terminal prior to running the command). The condor command is relatively general, but the slurm command assumes you are running on hipergator.
```
condor_submit_workers -M ${USER}-workqueue-coffea -t 900 --cores 4 --memory 18000 --disk 24000 1
slurm_submit_workers --cores 64 --memory 500000 -M ${USER}-workqueue-coffea -p "--partition hpg-default  --account avery --qos avery-b --time 2:00:00" 1
```

