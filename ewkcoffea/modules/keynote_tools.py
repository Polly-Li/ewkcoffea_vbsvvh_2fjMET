#!/bin/env python

import sys
import os
import math

###############################################################################################
#
#  Errors
#
###############################################################################################


# -*- coding: UTF-8 -*-

class E:
    """
    Properly propagates errors using all standard operations
    """
    def __init__(self, val, err=None):
        # assume poisson
        if err is None: err = abs(1.0*val)**0.5
        self.val, self.err = 1.0*val, 1.0*err

    def __add__(self, other):
        other_val, other_err = self.get_val(other)
        new_val = self.val + other_val
        new_err = (self.err**2.0 + other_err**2.0)**0.5
        return E(new_val, new_err)

    __radd__ = __add__

    def __sub__(self, other):
        other_val, other_err = self.get_val(other)
        new_val = self.val - other_val
        new_err = (self.err**2.0 + other_err**2.0)**0.5
        return E(new_val, new_err)

    def __rsub__(self, other):
        other_val, other_err = self.get_val(other)
        new_val = -(self.val - other_val)
        new_err = (self.err**2.0 + other_err**2.0)**0.5
        return E(new_val, new_err)


    def __mul__(self, other):
        other_val, other_err = self.get_val(other)
        new_val = self.val * other_val
        new_err = ((self.err * other_val)**2.0 + (other_err * self.val)**2.0)**0.5
        return E(new_val, new_err)

    __rmul__ = __mul__

    def __div__(self, other):
        other_val, other_err = self.get_val(other)
        new_val = self.val / other_val
        new_err = ((self.err/other_val)**2.0+(other_err*self.val/(other_val)**2.0)**2.0)**0.5
        return E(new_val, new_err)

    def __rdiv__(self, other):
        other_val, other_err = self.get_val(other)
        new_val = other_val / self.val
        new_err = ((other_err/self.val)**2.0+(self.err*other_val/(self.val)**2.0)**2.0)**0.5
        return E(new_val, new_err)

    def __pow__(self, other):
        # doesn't accept an argument of class E, only normal number
        new_val = self.val ** other
        new_err = ((other * self.val**(other-1) * self.err)**2.0)**0.5
        return E(new_val, new_err)

    def __neg__(self):
        return E(-1.*self.val, self.err)

    def __lt__(self, other):
        return self.val < other.val

    def get_val(self, other):
        other_val, other_err = other, 0.0
        if type(other)==type(self):
            other_val, other_err = other.val, other.err
        return other_val, other_err

    def round(self, ndec):
        if ndec == 0:
            self.val = int(self.val)
        else:
            self.val = round(self.val,ndec)
        self.err = round(self.err,ndec)
        return self

    def rep(self):
        use_ascii = False
        if use_ascii:
            sep = "+-"
        else:
            # sep = u"\u00B1".encode("utf-8")
            sep = u'\u00B1'
        if type(self.val).__name__ == "ndarray":
            import numpy as np
            # trick:
            # want to use numpy's smart formatting (truncating,...) of arrays
            # so we convert value,error into a complex number and format
            # that 1D array :)
            formatter = {"complex_kind": lambda x:"%5.2f {} %4.2f".format(sep) % (np.real(x),np.imag(x))}
            return np.array2string(self.val+self.err*1j,formatter=formatter, suppress_small=True, separator="   ")
        else:
            return "%s %s %s" % (str(self.val), sep, str(self.err))

    __str__ = rep

    __repr__ = rep

    def __getitem__(self, idx):
        if idx==0: return self.val
        elif idx==1: return self.err
        else: raise IndexError

r = None
def get_significance(exp,obs):
    """
    https://root.cern.ch/root/html526/RooStats__NumberCountingUtils.html
    """
    global r
    if not r: import ROOT as r
    return r.RooStats.NumberCountingUtils.BinomialObsZ(obs[0], exp[0], exp[1]/exp[0])

###############################################################################################
#
#  PyTable
#
###############################################################################################

def round_sig(x, sig=2):
    if x < 0.001: return x
    return round(x, sig-int(math.floor(math.log10(x)))-1)

class Table():

    def __init__(self):
        self.matrix = []
        self.colnames = []
        self.colsizes = []
        self.hlines = []
        self.rowcolors = {}
        self.extra_padding = 1
        self.d_style = {}
        self.set_theme_fancy()
        self.use_color = True

    def shorten_string(self, val, length):
        return val[:length//2-1] + "..." + val[-length//2+2:]
    def fmt_string(self, val, length, fill_char=" ", justify="c", bold=False, offcolor=False, color=None):
        ret = ""
        val = str(val)
        # lenval = len(val.decode("utf-8"))
        lenval = len(val)
        if lenval > length: val = self.shorten_string(val, length)
        if justify == "l":
            nr = (length-lenval-1)
            ret = " " + val + fill_char*nr
            # ret = " "+val.ljust(length-1, fill_char)
        elif justify == "r": ret = val.rjust(length, fill_char)
        elif justify == "c":
            nl = (length-lenval)//2
            nr = (length-lenval)-(length-lenval)//2
            ret = fill_char*nl + val + fill_char*nr
        if bold and self.use_color:
            ret = '\033[1m' + ret + '\033[0m'
        if offcolor and self.use_color:
            ret = '\033[2m' + ret + '\033[0m'
        if self.use_color:
            if color == "green":
                ret = '\033[00;32m' + ret + '\033[0m'
            if color == "blue":
                ret = '\033[00;34m' + ret + '\033[0m'
            if color == "lightblue":
                ret = '\033[38;5;117m' + ret + '\033[0m'
        return ret

    def set_theme_fancy(self):
        self.d_style["INNER_HORIZONTAL"] = '\033(0\x71\033(B'
        self.d_style["INNER_INTERSECT"] = '\033(0\x6e\033(B'
        self.d_style["INNER_VERTICAL"] = '\033(0\x78\033(B'
        self.d_style["OUTER_LEFT_INTERSECT"] = '\033(0\x74\033(B'
        self.d_style["OUTER_LEFT_VERTICAL"] = self.d_style["INNER_VERTICAL"]
        self.d_style["OUTER_RIGHT_INTERSECT"] = '\033(0\x75\033(B'
        self.d_style["OUTER_RIGHT_VERTICAL"] = self.d_style["INNER_VERTICAL"]
        self.d_style["OUTER_BOTTOM_HORIZONTAL"] = self.d_style["INNER_HORIZONTAL"]
        self.d_style["OUTER_BOTTOM_INTERSECT"] = '\033(0\x76\033(B'
        self.d_style["OUTER_BOTTOM_LEFT"] = '\033(0\x6d\033(B'
        self.d_style["OUTER_BOTTOM_RIGHT"] = '\033(0\x6a\033(B'
        self.d_style["OUTER_TOP_HORIZONTAL"] = self.d_style["INNER_HORIZONTAL"]
        self.d_style["OUTER_TOP_INTERSECT"] = '\033(0\x77\033(B'
        self.d_style["OUTER_TOP_LEFT"] = '\033(0\x6c\033(B'
        self.d_style["OUTER_TOP_RIGHT"] = '\033(0\x6b\033(B'

    def set_theme_basic(self):
        self.use_color = False

        self.d_style["INNER_HORIZONTAL"] = '-'
        self.d_style["INNER_INTERSECT"] = '+'
        self.d_style["INNER_VERTICAL"] = '|'
        self.d_style["OUTER_LEFT_INTERSECT"] = '|'
        self.d_style["OUTER_LEFT_VERTICAL"] = '|'
        self.d_style["OUTER_RIGHT_INTERSECT"] = '+'
        self.d_style["OUTER_RIGHT_VERTICAL"] = '|'
        self.d_style["OUTER_BOTTOM_HORIZONTAL"] = '-'
        self.d_style["OUTER_BOTTOM_INTERSECT"] = '+'
        self.d_style["OUTER_BOTTOM_LEFT"] = '+'
        self.d_style["OUTER_BOTTOM_RIGHT"] = '+'
        self.d_style["OUTER_TOP_HORIZONTAL"] = '-'
        self.d_style["OUTER_TOP_INTERSECT"] = '+'
        self.d_style["OUTER_TOP_LEFT"] = '+'
        self.d_style["OUTER_TOP_RIGHT"] = '+'

    def set_theme_csv(self):
        self.use_color = False

        self.d_style["INNER_HORIZONTAL"] = ''
        self.d_style["INNER_INTERSECT"] = ','
        self.d_style["INNER_VERTICAL"] = ','
        self.d_style["OUTER_LEFT_INTERSECT"] = ''
        self.d_style["OUTER_LEFT_VERTICAL"] = ''
        self.d_style["OUTER_RIGHT_INTERSECT"] = ''
        self.d_style["OUTER_RIGHT_VERTICAL"] = ''
        self.d_style["OUTER_BOTTOM_HORIZONTAL"] = ''
        self.d_style["OUTER_BOTTOM_INTERSECT"] = ''
        self.d_style["OUTER_BOTTOM_LEFT"] = ''
        self.d_style["OUTER_BOTTOM_RIGHT"] = ''
        self.d_style["OUTER_TOP_HORIZONTAL"] = ''
        self.d_style["OUTER_TOP_INTERSECT"] = ''
        self.d_style["OUTER_TOP_LEFT"] = ''
        self.d_style["OUTER_TOP_RIGHT"] = ''

    def set_theme_latex(self):
        self.d_style["INNER_HORIZONTAL"] = ''
        self.d_style["INNER_INTERSECT"] = ''
        self.d_style["INNER_VERTICAL"] = ' & '
        self.d_style["OUTER_LEFT_INTERSECT"] = ''
        self.d_style["OUTER_LEFT_VERTICAL"] = ''
        self.d_style["OUTER_RIGHT_INTERSECT"] = ''
        self.d_style["OUTER_RIGHT_VERTICAL"] = '\\\\ \\hline'
        self.d_style["OUTER_BOTTOM_HORIZONTAL"] = ''
        self.d_style["OUTER_BOTTOM_INTERSECT"] = ''
        self.d_style["OUTER_BOTTOM_LEFT"] = ''
        self.d_style["OUTER_BOTTOM_RIGHT"] = ''
        self.d_style["OUTER_TOP_HORIZONTAL"] = ''
        self.d_style["OUTER_TOP_INTERSECT"] = ''
        self.d_style["OUTER_TOP_LEFT"] = ''
        self.d_style["OUTER_TOP_RIGHT"] = ''


    def set_column_names(self, cnames):
        self.colnames = cnames
        self.update()

    def add_row(self, row, color=None):
        self.matrix.append(row)
        if color:
            self.rowcolors[len(self.matrix)] = color

    def add_column(self, colname, values):
        # if no matrix to begin with, just add the column, otherwise append to rows
        if len(self.matrix) == 0:
            for val in values:
                self.matrix.append([val])
        else:
            for irow in range(len(self.matrix)):
                if irow < len(values): val = values[irow]
                else: val = "-"
                self.matrix[irow].append(val)
        self.colnames.append(colname)


    def add_line(self):
        # draw hlines by making list of the row
        # indices, which we will check when drawing
        # the matrix
        self.hlines.append(len(self.matrix))

    def update(self):
        if not self.colnames:
            if self.matrix:
                self.colnames = range(1,len(self.matrix[0])+1)
        if self.matrix:
            for ic, cname in enumerate(self.colnames):
                self.colsizes.append(max(
                    max([len(str(r[ic])) for r in self.matrix])+2,
                    len(str(cname))+2))

    def sort(self, column=None, descending=True):
        self.update()
        icol = self.colnames.index(column)
        # sort matrix and range of numbers to get sorted indices for later use
        self.matrix, self.sortedidxs = zip( *sorted( zip( self.matrix,range(len(self.matrix))), key=lambda x: x[0][icol], reverse=descending))
        self.matrix = list(self.matrix)

        # now update row colors and hlines to match sorted matrix
        # one based indexing, not zero, so add 1
        oldtonewidx = dict([(self.sortedidxs[i]+1,i+1) for i in range(len(self.sortedidxs))])
        newrowcolors = {}
        for key,val in self.rowcolors.items():
            newrowcolors[oldtonewidx[key]] = val
        self.rowcolors = newrowcolors
        self.hlines = [oldtonewidx[hline] for hline in self.hlines]

    def print_table(self, **kwargs):
        print("".join(self.get_table_string(**kwargs)))

    def get_table_string(self, bold_title=True, show_row_separators=False, show_alternating=False, ljustall=False, show_colnames=True):
        self.update()
        nrows = len(self.matrix) + 1

        for irow,row in enumerate([self.colnames]+self.matrix):

            # line at very top
            if irow == 0:
                yield self.d_style["OUTER_TOP_LEFT"]
                for icol,col in enumerate(row):
                    yield self.d_style["OUTER_TOP_HORIZONTAL"]*(self.colsizes[icol]+self.extra_padding)
                    if icol != len(row)-1: yield self.d_style["OUTER_TOP_INTERSECT"]
                yield self.d_style["OUTER_TOP_RIGHT"]+"\n"

            if not show_colnames and irow == 0: continue

            # lines separating columns
            yield self.d_style["OUTER_LEFT_VERTICAL"]
            oc = False if not show_alternating else (irow%2==1 )
            bold = False if not bold_title else (irow==0)
            color = self.rowcolors.get(irow,None)
            if irow == 0: color = "lightblue"
            for icol,col in enumerate(row):
                j = "l" if icol == 0 else "c"
                if ljustall: j = "l"
                yield self.fmt_string(col, self.colsizes[icol]+self.extra_padding, justify=j, bold=bold,offcolor=oc,color=color)
                if icol != len(row)-1: yield self.d_style["INNER_VERTICAL"]
            yield self.d_style["OUTER_RIGHT_VERTICAL"]+"\n"

            # lines separating rows
            if (show_row_separators and (irow < nrows-1)) or (irow == 0):
                yield self.d_style["OUTER_LEFT_INTERSECT"]
                for icol,col in enumerate(row):
                    yield self.d_style["INNER_HORIZONTAL"]*(self.colsizes[icol]+self.extra_padding)
                    if icol != len(row)-1: yield self.d_style["INNER_INTERSECT"]
                yield self.d_style["OUTER_RIGHT_INTERSECT"]+"\n"

            # line at very bottom
            if irow == nrows-1:
                yield self.d_style["OUTER_BOTTOM_LEFT"]
                for icol,col in enumerate(row):
                    yield self.d_style["OUTER_BOTTOM_HORIZONTAL"]*(self.colsizes[icol]+self.extra_padding)
                    if icol != len(row)-1: yield self.d_style["OUTER_BOTTOM_INTERSECT"]
                yield self.d_style["OUTER_BOTTOM_RIGHT"]+"\n"
            else:
                # extra hlines
                if irow in self.hlines:
                    yield self.d_style["OUTER_LEFT_INTERSECT"]
                    for icol,col in enumerate(row):
                        yield self.d_style["OUTER_TOP_HORIZONTAL"]*(self.colsizes[icol]+self.extra_padding)
                        if icol != len(row)-1: yield self.d_style["INNER_INTERSECT"]
                    yield self.d_style["OUTER_RIGHT_INTERSECT"]+"\n"







#~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
#~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*~*
#~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=


#______________________________________________________________________________________________________________________
def human_format(num):
    is_fraction = False
    if num < 1:
        is_fraction = True
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    # add more suffixes if you need them
    if is_fraction:
        return '%.2g%s' % (num, ['', 'K', 'M', 'G', 'T', 'P'][magnitude])
    else:
        return '%.3g%s' % (num, ['', 'K', 'M', 'G', 'T', 'P'][magnitude])


#______________________________________________________________________________________________________________________
def yield_str(nominal, error, options={}):
    precuse = 4
    if "yield_prec" in options and type(options["yield_prec"]) is int:
        precuse = options["yield_prec"]
    noerror = False
    if "noerror" in options and options["noerror"]:
        noerror = options["noerror"]
    if precuse == 0 or noerror:
        return str(int(nominal))
    if noerror:
        return "{{:.{}g}}".format(precuse).format(nominal)
    else:
        e = E(nominal, error)
        if "human_format" in options:
            if options["human_format"]:
                # sep = u"\u00B1".encode("utf-8")
                sep = u'\u00B1'
                return "%s %s %s" % (human_format(e.val), sep, human_format(e.err))
            else:
                return e.round(precuse)
        else:
            # return e.round(precuse)
            # sep = u"\u00B1".encode("utf-8")
            sep = u'\u00B1'
            return "%s %s %s" % ('{{:.{}g}}'.format(precuse).format(e.val), sep, '{{:.{}g}}'.format(precuse).format(e.err))

#______________________________________________________________________________________________________________________
def check_regions_exists(yields_dict, reference, proc):
    if reference != yields_dict.keys():
        print("ewkcoffea.modules.keynote_tools.print_table(...): List of regions do not agree!")
        print("for process = {} there are following regions:", proc)
        for reg in yields_dict[proc].keys():
            print("  {}", reg)
        print("while for reference process there are following regions:")
        for reg in reference:
            print("  {}", reg)
        sys.exit()


#______________________________________________________________________________________________________________________
def get_keynote_str(table_str, dosimple=False):

    lines = table_str.split("\n")

    rtn_str = ""

    # delim=","
    delim=","

    if dosimple:
        for i, line in enumerate(lines):
            if i == 1:
                line = "".join([""+delim] + line.split()[3:])
                line = line.replace("|", "{}".format(delim))
                rtn_str += line[:-1] + "\n"
            elif i == 2:
                continue
            elif "|" in line:
                if len(line.split()) > 1:
                    title = line.split()[1] + delim
                else:
                    title = delim
                line = "".join([title] + line.split()[3:])
                line = line.replace("|", delim)
                # line = line.replace(u"\u00B1".encode("utf-8"), " " + u"\u00B1".encode("utf-8") + " ")
                line = line.replace(u"\u00B1", " " + u"\u00B1" + " ")
                rtn_str += line[:-1] + "\n"
    else:
        for i, line in enumerate(lines):
            if i == 1:
                line = "".join([""+delim] + line.split()[3:])
                line = line.replace("|", "{}{}{}".format(delim, delim, delim))
                rtn_str += line[:-1] + "\n"
            elif i == 2:
                continue
            elif "|" in line:
                if len(line.split()) > 1:
                    title = line.split()[1] + delim
                else:
                    title = delim
                line = "".join([title] + line.split()[3:])
                line = line.replace("|", delim)
                # line = line.replace(u"\u00B1".encode("utf-8"), delim + u"\u00B1".encode("utf-8") + delim)
                line = line.replace(u"\u00B1", delim + u"\u00B1" + delim)
                rtn_str += line[:-1] + "\n"

    return rtn_str

#______________________________________________________________________________________________________________________
def print_table(yields_dict, sigs=[], bkgs=[], yields_dict_data=None, data=None, systvar="nominal", options={}):
    # sanity checks
    if len(sigs) == 0 and len(bkgs) == 0 and data == None:
        print("ewkcoffea.modules.keynote_tools.print_table(...): You're not requesting anything to be printed! Check your arguments!")
        sys.exit()

    # if yields_dict_data is provided use that
    # if not use the same one as the other dict
    if yields_dict_data:
        yields_dict_data_ = yields_dict_data
    else:
        yields_dict_data_ = yields_dict

    # getting a list of regions and checking for all processes they all exists
    regions_ = list(yields_dict.keys())
    regions_.sort()
    # if data: check_regions_exists(yields_dict_data_, regions_, data)

    # restrict regions to pattern
    regions = []
    if "region_patterns" in options and options["region_patterns"]:
        # check that region_patterns is a list type
        if type(options["region_patterns"]) is not list:
            print("ewkcoffea.modules.keynote_tools.print_table(...): option 'region_patterns' is not a list type!")
            sys.exit()
        for i in regions_:
            for pattn in  options["region_patterns"]:
                if pattn in i:
                    regions.append(i)
                    break
    else:
        regions = regions_

    # If option "regions" is given provide the list specifically in that order
    if "regions" in options and options["regions"]:
        # check that regions is a list type
        if type(options["regions"]) is not list:
            print("ewkcoffea.modules.keynote_tools.print_table(...): option 'regions' is not a list type!")
            sys.exit()
        for r in options["regions"]:
            if r not in regions_:
                print(f"ewkcoffea.modules.keynote_tools.print_table(...): provided region={r} in 'regions' option is not present in the pkl file!")
                sys.exit()
        regions = options["regions"]

    x = Table()
    x.add_column("Regions", regions)

    for sig in sigs:
        x.add_column(sig, [ yield_str(yields_dict[reg][systvar][sig][0], math.sqrt(yields_dict[reg][systvar][sig][1]), options=options) for reg in regions ])

    for bkg in bkgs:
        x.add_column(bkg, [ yield_str(yields_dict[reg][systvar][bkg][0], math.sqrt(yields_dict[reg][systvar][bkg][1]), options=options) for reg in regions ])

    # total signal
    if len(sigs) != 0:
        totsig_yld_str = []
        for reg in regions:
            totb = E(0, 0)
            for sig in sigs:
                totb += E(yields_dict[reg][systvar][sig][0], math.sqrt(yields_dict[reg][systvar][sig][1]))
            totsig_yld_str.append(yield_str(totb.val, totb.err, options=options))
        x.add_column("Signal", totsig_yld_str)

    # total background
    if len(bkgs) != 0:
        totbkg_yld_str = []
        totbs = {}
        for reg in regions:
            totb = E(0, 0)
            for bkg in bkgs:
                totb += E(yields_dict[reg][systvar][bkg][0], math.sqrt(yields_dict[reg][systvar][bkg][1]))
            totbs[reg] = totb
            totbkg_yld_str.append(yield_str(totb.val, totb.err, options=options))
        x.add_column("Background", totbkg_yld_str)

    # Add data columns
    if data:
        x.add_column("Data", [ yield_str(yields_dict_data_[reg][systvar][data][0], math.sqrt(yields_dict_data_[reg][systvar][data][1]), options=options) for reg in regions ])

        # Also add data/MC
        data_mc_str = []
        data_mcs = {}
        dyields = {}
        for reg in regions:
            dyield = E(yields_dict_data_[reg][systvar][data][0], math.sqrt(yields_dict_data_[reg][systvar][data][1]))
            dyields[reg] = dyield
            data_mcs[reg] = dyields[reg].__div__(totbs[reg])
            data_mc_str.append(yield_str(data_mcs[reg].val, data_mcs[reg].err, options=options))

        x.add_column("Data/MC", data_mc_str)

    x.print_table()
    x.set_theme_basic()

    output_name = "yields"

    if "output_name" in options and options["output_name"]:
        if type(options["output_name"]) is not str:
            print("ewkcoffea.modules.keynote_tools.print_table(...): option 'output_name' is not a str type!")
            sys.exit()
        output_name = options["output_name"]

    # Write text version
    os.makedirs("tables", exist_ok=True)
    f = open(f"tables/{output_name}.txt", "w")
    f.write("".join(x.get_table_string()))

    f = open(f"tables/{output_name}.csv", "w")
    if "yield_prec" in options and type(options["yield_prec"]) is int:
        f.write(get_keynote_str("".join(x.get_table_string()), dosimple=True if options["yield_prec"] == 0 else False))
    else:
        f.write(get_keynote_str("".join(x.get_table_string())))

