import pysam
from collections import Counter
import argparse
from os.path import isfile
from itertools import chain

record = []
c = Counter()
endResult = []

with open("merge_pcr.txt", "w"):
     pass

def flag_check(bam_file):
    record = []
    sam = pysam.AlignmentFile(bam_file, "rb")
    iter = sam.fetch(multiple_iterators=True, until_eof=True)
    for line in iter:
        flag_seq = []
        flag = line.flag
        tags = list(line.get_tags())
        print(line)
        if (line.is_unmapped == False and 'XS' not in chain(*tags)):
            if flag == 0:
                flag_seq.append('+')
            elif line.is_reverse:
                flag_seq.append('-')
            record.append((line.reference_start,
                           line.reference_start + line.reference_length,
                           line.query_name,
                           flag_seq,
                           sam.get_reference_name(line.reference_id),
                           line.get_reference_sequence))
    sam.close()
    record = sorted(record, key=lambda x: x[2])
    return record

def fastq_read(fastq_file):
    fastaq = []
    with open(fastq_file, 'r') as file:
        r = file.readlines()
        for y in r:
            if r.index(y) % 4 == 0:
                b = (str(y)).strip().split()
                bar_code = (r[(r.index(y)) + 1].strip())
                fastaq.append((b[0][1:], bar_code))

    fastaq = sorted(fastaq, key=lambda x: x[0])
    return fastaq

rec = []
pos = []
def chromosome_count(i,j):
    s_range = (record[i][0], record[i][1])
    c.update([s_range])
    if (s_range not in pos):
        rec.append((record[i][4], record[i][0], record[i][1], record[i][2], ",".join((record[i][3]))))
        pos.append(s_range)
    if len(fastaq) - 1 == j:
        pass
    else:
        j += 1
    return rec, j

def chromosome_counter():
    j =0
    reco = []

    for i in range(len(record)):
        if record[i][2] == fastaq[j][0]:
            reco, j = chromosome_count(i, j)
        elif record[i][2] > fastaq[j][0]:
            while record[i][2] > fastaq[j][0]:
                j += 1
            if record[i][2] == fastaq[j][0]:
                reco, j = chromosome_count(i, j)
        else:
            while record[i][2] < fastaq[j][0]:
                continue
    return reco

def printing(endResult):
    with open("merge_pcr.txt", "w") as f:
        for entry in endResult:
            wr = ("%s\t %s\t %s\t %s\t %s\t %s" % (entry[0], entry[1], entry[2], entry[3], c[(entry[1],entry[2])], entry[4]))
            f.write(wr)
            f.write('\n')
            print("%s\t %s\t %s\t %s\t %s\t %s" % (entry[0], entry[1], entry[2], entry[3], c[(entry[1],entry[2])], entry[4]))

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
parser.add_argument("bam_file", help="Path to bam file containing alignments.", metavar='BAM File')
parser.add_argument("fastq_file", help="Path to fastq barcode library.", metavar='FASTQ File')
parser.add_argument("-o", "--outfile", required=False, help="Write results to this file.",
                    metavar='Output File')
args = parser.parse_args()
record = flag_check(args.bam_file)
fastaq = fastq_read(args.fastq_file)

try:
    endResult = chromosome_counter()


except:
    if not isfile(args.fastq_file):
        print("ERROR: Fastq'{}' not found.".format(args.fastq_file))
    if not isfile(args.bam_file):
        print("ERROR: bam file '{}' not found.".format(args.bam_file))
    if not isfile(args.out_file):
        print("ERROR: bam file '{}' not found.".format(args.out_file))
    exit()
printing(endResult)
