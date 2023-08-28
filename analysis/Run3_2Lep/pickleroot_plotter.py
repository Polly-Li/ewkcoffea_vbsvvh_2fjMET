import argparse
import gzip
import sys
import pickle
import ROOT as rt

parser = argparse.ArgumentParser(description='You can customize your run')
parser.add_argument('--prefix', '-r'   , nargs='?', default='', help = 'Prefix to look for the files')
parser.add_argument('pickleFiles'        , nargs='+', default='', help = 'Input file(s) containing hists')
parser.add_argument('--outname','-o'   , default='plotsTopEFT', help = 'Name of the output file with histograms')
parser.add_argument('--do_DataMerge'      , action='store_true', help = 'Merge Data plots (Default is False)')


args = parser.parse_args()
prefix = args.prefix
pickleFiles  = args.pickleFiles
outname    = args.outname
do_DataMerge  = args.do_DataMerge

#Do some basic input checks
print("Input Filie(s): ",pickleFiles)
if prefix != '':
    print("All will accessed with prefix ",prefix)
if (len(pickleFiles) < 2) and do_DataMerge:
    print("ERROR: You have asked to merge hists, but have only given one file!")
    sys.exit()
if (len(pickleFiles) > 1) and (not do_DataMerge):
    print("ERROR: You have given more than one file for inputs, but did not ask to merge hists!")
    sys.exit()
if (len(pickleFiles) < 1):
    print(ERROR: No Input Files given!)
    sys.exit()


#Create output
o = rt.TFile(f"{outname}.root", "recreate");
o.cd()


#Create the file path(s)
file_path = [prefix + item for item in pickleFiles]

#Seperate the cases for multiple input files

#-------------------1 Input File -----------------------------------------------------------------
if len(pickleFiles) == 1:
    with gzip.open(file_path[0],'rb') as f:
        data = pickle.load(f)
    
    hist_lst = list(data.keys())
    
    for hist_name in hist_lst:
        pickle_hist = data[hist_name]
        bin_edges = list(pickle_hist.axes[0].edges)
        nbins = pickle_hist.axes[0].size
        root_hist = rt.TH1F(f"{hist_name}",f"{hist_name}",nbins,bin_edges[0],bin_edges[-1])
        
        for i in range(1, nbins+1): #ROOT hist binning starts at 1 
            root_hist.SetBinContent(i,pickle_hist[i-1])
        root_hist.Write()    
#            
#
#print("File has been created. Happy Plotting!")
