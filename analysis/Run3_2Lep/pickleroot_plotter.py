#Stuff here
import argparse
import gzip
import sys
import pickle
import ROOT as rt

parser = argparse.ArgumentParser(description='You can customize your run')
parser.add_argument('--prefix', '-r'     , nargs='?', default='', help = 'Prefix to look for the files')
parser.add_argument('pickleFiles'        , nargs='+', default='', help = 'Input file(s) containing hists')
parser.add_argument('--outname','-o'     , nargs='*', default=[], help = 'Output file(s) names (Deafult is the same name as input)')
#parser.add_argument('--outname','-o'     , default='default_name', help = 'Name of the output file with histograms')


args = parser.parse_args()
prefix = args.prefix
pickleFiles  = args.pickleFiles
print(pickleFiles)
outname = args.outname
if len(outname) == 0:
    outname = [item.replace('.pkl.gz', '') for item in pickleFiles]
print(outname)


#Do some basic input checks
print("Input Filie(s): ",pickleFiles)
if prefix != '':
    print("All will accessed with prefix ",prefix)
if (len(pickleFiles) < 1):
    print("ERROR: No Input Files given!")
    sys.exit()
if len(pickleFiles) != len(outname):
    print("ERROR: Input and Output File lengths do not match!")


#Create the file path(s)
file_path = [prefix + item for item in pickleFiles]

for i in range(len(pickleFiles)):
    #Create output
    o = rt.TFile(f"{outname[i]}.root", "recreate");
    o.cd()

    with gzip.open(file_path[i],'rb') as f:
        data = pickle.load(f)
    
    hist_lst = list(data.keys())
    
    for hist_name in hist_lst:
        pickle_hist = data[hist_name]
        bin_edges = list(pickle_hist.axes[0].edges)
        nbins = pickle_hist.axes[0].size
        root_hist = rt.TH1F(f"{hist_name}",f"{hist_name}",nbins,bin_edges[0],bin_edges[-1])
        
        for j in range(1, nbins+1): #ROOT hist binning starts at 1
            root_hist.SetBinContent(j,pickle_hist[j-1].value)
        root_hist.Write()    
            

print("File has been created. Happy Plotting!")
