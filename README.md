# ewkcoffea
This analysis repository contains scripts and tools for performing analyses associated with EWK physics within the `coffea` framework, making use of the `topcoffea` package. 

## Setup instructions

First, clone the repository and `cd` into the toplevel directory. 
```
https://github.com/kmohrman/ewkcoffea.git
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
