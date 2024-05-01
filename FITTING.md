# Statistical analysis

This readme describes how to perform the statistical analysis on the histograms produced by `ewkcoffea`.

## Making datacards 

The information in the histograms (produced by the processor) can be converted into datacards (to be fed to `combine`) via the `make_datacards.py` script.

* Example of running over Run2 histograms, with data-driven background estimation (on by default) and systematics (`-s`) for the BDT bins (`--bdt`): 
  ```
  python make_datacards.py histos/yourhist.pkl.gz --run run2 -s --bdt
  ```

* Example of running over Run3 (which right now is just 2022) histograms: 
  ```
  python make_datacards.py histos/yourhisto.pkl.gz --run run3
  ```

## Using combine to perform the statistical analysis 

The statistical analysis is performed with the `combine` tool. This should be set up _outside_ of your conda env, following the instructions [here](https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/latest/#within-cmssw-recommended-for-cms-users). 

Copy the datacards (from `ewkcoffea`) to your `combine` working area. To combine the cards from all of the channels into a single card, use the `combineCards.py` command, e.g. run the following (in a directory where `wwz4l_card_*` match to all of your cards and no extra cards): 
```
combineCards.py wwz4l_card_* > comb.txt
```
To check the significance, run the following:
```
combine -M Significance comb.txt --expectSignal=1 -t -1
```
