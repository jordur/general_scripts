#!/usr/bin/env python

import csv
import argparse
import collections


def change_type(l, dtype=float):
    return map(dtype, l)

parser = argparse.ArgumentParser(description='Normalize consam CSV by the number of reads.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('csv', type=str,
                    help='input: CSV to read chr,value1,value2.')

parser.add_argument('normalized_csv', type=str,
                    help='output: CSV normalized by number reads.')

args = parser.parse_args()

# Prepare the list of chromosomes
chromosomes = collections.OrderedDict()
for chromosome in range(1, 23):
    chromosomes[str(chromosome)] = []
chromosomes['X'] = []
chromosomes['Y'] = []

total_reads = 0

with open(args.csv, 'rb') as f:
    reader = csv.reader(f)
    for row in reader:
        if row[0] not in ('X', 'Y'):
            float_row = [float(i) for i in row[1:]]
            total_reads += sum(float_row)

    f.seek(0)

    for row in reader:
        if row[0] not in ('X', 'Y'):
            float_row = [float(i) for i in row[1:]]
            normalized_row = [value/total_reads for value in float_row]
            chromosomes[row[0]].extend(normalized_row[0:])
        else:
            float_row = [float(i) for i in row[1:]]
            chromosomes[row[0]].extend(float_row[0:])

writer = csv.writer(open(args.normalized_csv, 'wb'), delimiter=',')
#for key in sorted(chromosomes.iterkeys()):
for key, value in chromosomes.items():
    #print key
    results_row = [key]
    results_row.extend(chromosomes[key])
    #string_values = [str(i) for i in value[0:]]
    #string_values = ','.join(str(elt) for elt in string_values)
    writer.writerow(results_row)