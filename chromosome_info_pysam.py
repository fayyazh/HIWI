import pysam
from collections import Counter
import argparse
from os.path import isfile
from Bio import SeqIO


bam_data = []
score = Counter()
endResult = []
merge_data = []
merge_index = []
fastq_data = {}

def bam_reader(bam_file):
    bam_filter = []
    bam = pysam.AlignmentFile(bam_file, "rb")
    data = bam.fetch(multiple_iterators=True, until_eof=True)

    for line in data:
        if (line.is_unmapped == False and line.has_tag("XS") == False):
            if line.is_reverse:
                flag = '-'
            else:
                flag = '+'
            bam_filter.append((line.reference_start,
                           line.reference_start + line.reference_length,
                           line.query_name, flag,
                           bam.get_reference_name(line.reference_id)))
    bam.close()
    return bam_filter

def fastq_reader(fastq_file):
    fastq_dt = {}
    for fastq_data in SeqIO.parse(fastq_file, "fastq"):
        fastq_dt.update({fastq_data.id : fastq_data.seq})
    return fastq_dt

def chromosome_counter(i):
    # merge_checks contains reference chromosome, ref. start, ref. end and strand.
    merge_checks = (bam_data[i][4], bam_data[i][0], bam_data[i][1], fastq_data[bam_data[i][2]], "".join(bam_data[i][3]))
    score.update([merge_checks])
    if(merge_checks not in merge_index):
        merge_data.append((bam_data[i][4], bam_data[i][0], bam_data[i][1], bam_data[i][2], ",".join((bam_data[i][3]))))
        merge_index.append(merge_checks)
    return merge_data

def chromosome_info():
    chr_info = []
    for i in range(len(bam_data)):
        rec_id = bam_data[i][2]
        if rec_id in fastq_data:
            chr_info = chromosome_counter(i)
        else:
            print(rec_id + " ID not found in fastq file")
    return chr_info

def printing(endResult):
    with open(args.output_file, "w") as f:
        for entry in endResult:
            wr = ("%s\t%s\t%s\t%s\t%s\t%s" % (entry[0], entry[1], entry[2], entry[3], score[(entry[0], entry[1], entry[2], fastq_data[entry[3]], entry[4])], entry[4]))
            f.write(wr + '\n')

tool_description = """
Merge PCR duplicates according to random barcode library.
Barcodes containing uncalled base 'N' are removed.
Input:
* bam file containing fastq read-id and other details
* fastq library of random barcodes
Output:
* bed6 file with random barcode in name field and number of PCR duplicates as
  score, sorted by fields chrom, start, stop, strand, name
Example usage:
python3 pysam_test.py example1.bam example2.fa
"""

epilog = """
Author: Fayyaz Hussain
Status: Testing
"""

# parse command line arguments
parser = argparse.ArgumentParser(prog="Chromosomes Information", description=tool_description,
                                 epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("bam_file", help="Path to bam file containing alignments.", metavar='BAM_File')
parser.add_argument("fastq_file", help="Path to fastq barcode library.", metavar='FASTQ_File')
parser.add_argument("-o", "--output_file", required=True, help="Write results to this file.",
                    metavar='Output_File')
args = parser.parse_args()
fastq_data = fastq_reader(args.fastq_file)
bam_data = bam_reader(args.bam_file)

try:
    endResult = chromosome_info()
except:
    if not isfile(args.fastq_file):
        print("ERROR: Fastq'{}' not found.".format(args.fastq_file))
    if not isfile(args.bam_file):
        print("ERROR: bam file '{}' not found.".format(args.bam_file))
    if not isfile(args.outfile):
        print("ERROR: bam file '{}' not found.".format(args.output_file))
    exit()
printing(endResult)
