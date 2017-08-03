#!/usr/bin/env python

import csv
import pickle
import argparse


def change_type(l, dtype=float):
    return map(dtype, l)

parser = argparse.ArgumentParser(description='Translate a csv to a WISECONDOR pickle file',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('csv', type=str,
                    help='input: csv to read chr,value1,value2,...')

parser.add_argument('pickle', type=str,
                    help='output: pickle to dump file')

args = parser.parse_args()

# Prepare the list of chromosomes
chromosomes = dict()
for chromosome in range(1, 23):
    chromosomes[str(chromosome)] = []
chromosomes['X'] = []
chromosomes['Y'] = []

with open(args.csv, 'rb') as f:
    reader = csv.reader(f)
    for row in reader:
        float_row = [float(i) for i in row[1:]]
        chromosomes[row[0]].extend(float_row[0:])

pickle.dump(chromosomes, open(args.pickle, 'wb'))

