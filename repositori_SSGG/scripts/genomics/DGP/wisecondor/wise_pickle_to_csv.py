#!/usr/bin/env python

import csv
import pickle
import argparse

parser = argparse.ArgumentParser(description='Translate a pickle to a csv file',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('pickle', type=str,
                    help='input: pickle file')
parser.add_argument('csv', type=str,
                    help='output: csv to dump table into')
args = parser.parse_args()

sample = pickle.load(open(args.pickle, 'rb'))
outWriter = csv.writer(open(args.csv, 'wb'))

keys = [str(x) for x in range(1, 23)]
keys.extend(['X', 'Y'])
for key in keys:
    row = [key]
    row.extend(sample[key])
    outWriter.writerow(row)