import uproot
import argparse
import gzip
import sys
import pickle

parser = argparse.ArgumentParser(description='You can customize your run')
parser.add_argument('--prefix', '-r'     , nargs='?', default='', help = 'Prefix to look for the files')
parser.add_argument('pickleFiles'        , nargs='+', default='', help = 'Input file(s) containing hists')
parser.add_argument('--outname','-o'     , nargs='*', default=[], help = 'Output file(s) names (Deafult is the same name as input)')

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
    with gzip.open(file_path[i],'rb') as f:
        data = pickle.load(f)

    hist_lst = list(data.keys())

    with uproot.recreate(f"{outname[i]}.root") as root_file:
        for hist_name in hist_lst:
            pickle_hist = data[hist_name]
            root_file[f"{hist_name}"] = pickle_hist
