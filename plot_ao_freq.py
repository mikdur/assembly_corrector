#!/usr/bin/env python3

import sys, argparse

# For reading vcf files
import vcf 

from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

parser = argparse.ArgumentParser()

parser.add_argument("--vcf", required=True, nargs="?", help="VCF file to be plotted")
parser.add_argument("--out", nargs="?", required=True,
                    help="Output PDF")
parser.add_argument("--title", nargs="?", default="",
                    help="Plot title")
parser.add_argument("--dp", nargs="?", type=float, default=0,
                    help="Minimum DP of a variant to use it for correction")
args=parser.parse_args()



# Read VCF file
vcf_reader = vcf.Reader(open(args.vcf, "r"))

# Extract relevant vcf data
poly_data = [ [r.INFO['DP'], r.INFO['AO'][0]] for r in vcf_reader if r.INFO['DP'] >= args.dp ]

# Allele frequency
a_f = [ float(x[1]) / float(x[0]) for x in poly_data ]


# Do the plotting

# Open output file
pp = PdfPages(args.out)

#plt.yscale('log')
plt.hist(a_f, bins=40, normed=False)
plt.xlabel("Alternate allele proportion")
plt.ylabel("Frequency")
plt.title(( args.vcf if args.title == "" else args.title) + " (n=%d)" % (len(a_f)))

plt.show()
pp.savefig()
pp.close()


