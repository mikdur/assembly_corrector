#!/usr/bin/env python3

import vcf, argparse, sys
from Bio import SeqIO, Seq


parser = argparse.ArgumentParser()

parser.add_argument("--vcf", required=True, nargs="?", help="VCF file with variations that should be used " +
                   "to correct the genome sequence", type=argparse.FileType('r'))
parser.add_argument("--min-dp", nargs="?", type=float, default=0,
                    help="Minimum DP of a variant to use it for correction")
parser.add_argument("--min-qual", nargs="?", type=float, default=0,
                    help="Minium quality of variant to use it for correction")
parser.add_argument("--min-ao-freq", nargs="?", type=float, default=0,
                    help="Minimum relative frequency of alternate observation")
parser.add_argument("--min-paired", nargs="?", type=float, default=0,
                    help="Minimun proportion of alternate allele supported by proper pairs")
parser.add_argument("--fasta-in", nargs="?", required=True,
                    help="Input fasta")
parser.add_argument("--fasta-out", nargs="?", required=True,
                    help="Output fasta")
parser.add_argument("--verbose", default=False, help="Be verbose", action='store_true')
args=parser.parse_args()


debug_mode = False
def debug(msg):
    if debug_mode:
        print(msg, file=sys.stderr)

if args.verbose:
    debug_mode = True

debug("Verbose mode enabled")


debug("Opening VCF file")
vcf_reader = vcf.Reader(args.vcf)

corrections = { }
for record in vcf_reader:
    # Get the alternate allele with the highest score
    ind = 0 if len(record.ALT) == 1 \
      else (sorted(zip(record.INFO['QA'], range(len(record.ALT))),
                 key=lambda a: a[0])[0][1])

    if record.INFO['DP'] > args.min_dp and \
      record.QUAL >= args.min_qual and \
      record.INFO['AO'][ind] / record.INFO['DP'] >= args.min_ao_freq and \
      record.INFO['PAIRED'][ind] >= args.min_paired:
        c = { "pos" : record.POS,
              "ref" : record.REF,
              "alt" : str(record.ALT[ind]) }
      
        if record.CHROM not in corrections:
            corrections[record.CHROM] = [ c ]
        else:
            corrections[record.CHROM].append(c)
            
print("Found %d corrections spread across %d scaffolds" % \
      (sum( [ len(v) for k,v in corrections.items() ] ),
       len(corrections)))

corrections = { k : sorted(v, key=lambda k: k['pos'], reverse=True) for k, v in corrections.items() }


print("Correcting...")
seq_records = SeqIO.parse(args.fasta_in, "fasta")
out_seqs = [ ]
for rec in seq_records:
    if rec.id in corrections:
        seq = str(rec.seq)

        debug("Correcting %s (%d bp, %d corrections)" % (rec.id, len(seq), len(corrections[rec.id])))

        for corr in corrections[rec.id]:

            if seq[corr['pos'] - 1 :corr['pos'] + len(corr['ref']) - 1].upper() != corr['ref'].upper():
                raise Exception("Reference of sequence file and VCF differs!")
            seq = seq[:corr['pos'] - 1] + corr['alt'] + seq[corr['pos'] + len(corr['ref']) - 1:]
        rec.seq = Seq.Seq(seq)
        out_seqs.append(rec)
    else:
        out_seqs.append(rec)

SeqIO.write(out_seqs, args.fasta_out, "fasta")
