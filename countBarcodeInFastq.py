
import io
import os
import sys
import gzip
import argparse
import numpy as np
import pandas as pd
import multiprocessing as mp
from functools import partial

def main():
    parser = argparse.ArgumentParser(
        description='Scanning a gzipped fastq to count barcodes '\
        'provided in an Excel file' )
    parser.add_argument('fastq', help='input fastq.gz')
    parser.add_argument('result', help='result file')
    parser.add_argument('-b', dest='barcode_setting', required=True,
                        help='barcode setting in excel')
    parser.add_argument('-p', dest='processes', help='processes number',
                        type=int)
    parser.add_argument('-s', dest='insert_size', help='library insert size',
                        type=int, required=True)
    parser.add_argument('-d', dest='direction', type=int,
                        required=True, help='library direction: '\
                        '{0 forward 1 backward 2 both}')
    args = parser.parse_args()
    barcode_combination = read_barcode_list_excel(args.barcode_setting,
                                                  args.insert_size)
    result = count_barcode_in_fastq(args.fastq, barcode_combination,
                                    identify_barcode_in_read, args.direction,
                                    args.processes)
    result = barcode_count_result_rearrange(result)
    print("Reads scanned: {}".format(result["total_reads"]))
    print("Reads with barcodes: {}".format(result["reads_identified"]))
    df = result['barcode_counts']
    print(df)
    writer = pd.ExcelWriter(args.result)
    df.to_excel(writer, index=False)
    writer.save()


def read_barcode_list_excel(barcode_excel, insert_size=100):
    barcode_dfs = pd.read_excel(barcode_excel, None)
    df_barcode_setting = barcode_dfs['barcode_position']
    barcode_sets = []
    for idx, row in df_barcode_setting.iterrows():
        barcode_set_name = row['barcode_set']
        if barcode_set_name not in barcode_dfs.keys():
            raise ValueError("No sheet '{}'".format(barcode_set_name))
        bs = None
        bs = BarcodeSet(barcode_set_name, int(row["start"]),
                        int(row["end"]), row["barcode_type"],                                           insert_size)
        for bcIdx, barcode_row in barcode_dfs[barcode_set_name].iterrows(): 
            bs.add(barcode_row['barcode_id'], barcode_row['barcode_sequence'])
        barcode_sets.append(bs)
    barcode_sets.sort(key=lambda x: x.start_position())
    return BarcodeCombination(barcode_sets)

def count_barcode_in_fastq(fastq_gz, barcode_combination, idenfunc,
                           direction=0, processes=1):

    reads, results, barcode_counts = [], [], []
    read_count = 0
    # attach arguments to barcode identification function 
    # (barcode_combination, library direction)
    iden = partial(idenfunc, barcode_combination=barcode_combination,
                   direction=direction)
    def count_barcode_in_reads(reads, processes, iden):
        if processes > 1:
            pool = mp.Pool(processes=processes)
            result = pool.map(iden, reads)
            pool.close()
            pool.join()
        else:
            result = [iden(read) for read in reads]
        return result

    with io.BufferedReader(gzip.open(fastq_gz, 'rb')) as f:
        print("\n" + "Counting barcodes in "\
              "{}...".format(os.path.basename(fastq_gz))) 
        line_count, read_count = 0, 0
        for line in f:
            line_count +=1
            if line_count % 4 != 2:  # skip non sequence line
                continue
            reads.append(line.decode()[:-1]) # add sequence into reads list
            read_count += 1
            if read_count % 100000 == 0:
                result = count_barcode_in_reads(reads, processes, iden)
                results.extend(result)
                barcode_counts.append(np.apply_along_axis(sum, 0,                                                             np.array(results)))
                reads = []
                results= []
                print('.', end='', flush=True)

            if len(barcode_counts) == 5:
                barcode_counts = [np.apply_along_axis(sum, 0, 
                          np.array(barcode_counts))]
    if reads:
        result = count_barcode_in_reads(reads, processes, iden)
        results.extend(result)
        barcode_counts.append(np.apply_along_axis(sum, 0, np.array(results)))
        print('.', end='', flush=True)
    print()
    barcode_counts = np.apply_along_axis(sum, 0, np.array(barcode_counts))
    return {'total_reads': read_count,
            'barcode_counts': barcode_counts, 
            'barcode_types': barcode_combination.barcode_combination_types(),
            'barcode_labels': barcode_combination.barcode_combination_labels()}

def barcode_count_result_rearrange(barcode_count_result):
    df = pd.DataFrame([barcode + [count] for barcode, count in 
             zip(barcode_count_result['barcode_labels'],
                 barcode_count_result['barcode_counts'])])
    for key, val in barcode_count_result['barcode_types'].items():
        df[key] = df.iloc[:,val[0]]
        for barcode_order, pos in enumerate(val):
            if len(val) > 1:
                df[key] = df[key] + '_' + df.iloc[:,pos]
    df['count'] = df.iloc[:,len(barcode_count_result['barcode_labels'][0])]
    return {'total_reads': barcode_count_result['total_reads'],
            'barcode_counts': df[list(barcode_count_result['barcode_types'].keys()) + ['count']],
            'reads_identified': df['count'].sum()
            }

def identify_barcode_in_read(read, barcode_combination, direction=0):
    """Identify occurence of barcode in a DNA sequence (read)
        Args:
            read: A DNA sequence consist of 'A', 'C', 'G', 'T' and 'N'
            barcode_combination: An BarcodeCombination object
        Returns:
            A 0,1 numpy array showing occurrence in read for each
            barcode combination 
    """

    # Get barcode information from BarcodeCombination Object
    barcode_positions = barcode_combination.barcode_combination_positions()
    barcode_combinations = barcode_combination.barcode_combination_in_number() 

    # Retrieve sequence from read at barcode positions
    sequences_at_barcode_positions = [ pattern_to_number(read[pos[0]:pos[1]])
                                      for pos in barcode_positions ]
    # Check occurrence with np.array_qual and np.apply_along_axis
    occurrence = np.apply_along_axis(
          lambda x: np.array_equal(x, 
                                    np.array(sequences_at_barcode_positions)),
                 1, barcode_combinations) + 0
    if direction == 0:  
        return occurrence
    
    # Get reverse barcode information from BarcodeCombination Object
    reverse_barcode_positions = barcode_combination.reverse_barcode_combination_positions()
    reverse_barcode_combinations = barcode_combination.reverse_barcode_combination_in_number()

    # Identify occurrence of barcode in read at reverse direction
    sequences_at_reverse_barcode_positions = [pattern_to_number(
          read[pos[0]:pos[1]]) for pos in
          reverse_barcode_positions  ]
    reverse_occurrence = np.apply_along_axis(
          lambda x: np.array_equal(x, 
                 np.array(sequences_at_reverse_barcode_positions)), 1, 
                 reverse_barcode_combinations) +0
    if direction == 1:
        return reverse_occurrence
    elif direction == 2:
        return occurrence + reverse_occurrence

def combinations(barcode_lists):
    if len(barcode_lists) == 1:
        if isinstance(barcode_lists[0][0], list):
            return barcode_lists[0]
        else:
            return [[x] for x in barcode_lists[0]]
    prefix_list = barcode_lists[0]
    if not isinstance(prefix_list[0], list):
        prefix_list = [[x] for x in barcode_lists[0]]
    suffix_list = barcode_lists[1]
    target_list = []
    for prefix in prefix_list:
        for suffix in suffix_list:
            target_list.append(prefix[:] + [suffix])
    if len(barcode_lists) == 2:
        return target_list
    else:
        return combinations([target_list]+barcode_lists[2:])

def pattern_to_number(pattern):
    symbol_dict = {'A':0, 'C':1, 'G':2, 'T':3, 'N':4}
    if len(pattern) == 1:
        return symbol_dict[pattern]
    return symbol_dict[pattern[0]] * (5**(len(pattern)-1)) + \
                    pattern_to_number(pattern[1:])

def reverse_complement(pattern):
    symbol_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
    return ''.join([symbol_dict[base] for base in pattern][::-1])



class DuplicateError(Exception):
    pass
    '''def __init__(self, message, errors):
        super().__init__(message)
        self.errors = errors
    '''


class BarcodeCombination:

    def __init__(self, barcode_sets):
        self._barcode_combination_types = [x.barcode_type() for x in
                                           barcode_sets]
        self._barcode_combination_types_dict = self.list2dict(
              self._barcode_combination_types)
        self._barcode_combination_labels = combinations( 
              [x.barcode_labels() for x in barcode_sets])
        self._barcode_combination_sequences = combinations(
              [x.barcode_sequences() for x in barcode_sets])
        self._barcode_combination_in_number = np.array(combinations(
              [x.barcode_in_number() for x in barcode_sets]))
        self._barcode_combination_positions = [(x.start_position(),
                                                x.end_position()) for x in
                                               barcode_sets]
        self._reverse_barcode_combination_positions =[
              (x.reverse_start_position(), x.reverse_end_position()) for x in
              barcode_sets]
        self._reverse_barcode_combination_in_number = np.array(combinations(
              [x.reverse_barcode_in_number() for x in barcode_sets ]))

    def list2dict(self, l):
        my_dict = {}
        for i, item in enumerate(l):
            my_dict[item] = my_dict.get(item, []) + [i]
        return my_dict

    def barcode_combination_positions(self):
        return self._barcode_combination_positions
    
    def reverse_barcode_combination_positions(self):
        return self._reverse_barcode_combination_positions

    def barcode_combination_sequences(self):
        return self._barcode_combination_sequences

    def barcode_combination_in_number(self):
        return self._barcode_combination_in_number

    def reverse_barcode_combination_in_number(self):
        return self._reverse_barcode_combination_in_number

    def barcode_combination_labels(self):
        return self._barcode_combination_labels

    def barcode_combination_types(self):
        return self._barcode_combination_types_dict



class BarcodeSet:
    """Return a barcode Set 
       
       Args:
           barcode_set_name: string
           start_position: start position of barcode in library
           end_position: end position of barcode in library
           insert_size: length of library to be sequenced
      
       Return:
           A barcode set with multiple attributes:
               barcode_set_name
               start_position
               end_position
               barcode_labels: A list
               barcode_sequences: a list
               barcode_in_number: conversion barcode sequences into integers
               reverse_start_position
               reverse_end_position
               reverse_barcode_in_number
    """

    def __init__(self, barcode_set_name,start_position, end_position,
                 barcode_type, insert_size):
        self._barcode_set_name = barcode_set_name
        self._barcode_type = barcode_type
        self._insert_size = insert_size
        self._start_position = start_position - 1
        self._end_position = end_position
        self._barcode_length = self._end_position - self._start_position
        self._reverse_start_position = insert_size - end_position
        self._reverse_end_position = insert_size - (start_position-1)
        self._barcode_sequences = []
        self._barcode_labels = []
        self._barcode_in_number = []
        self._reverse_barcode_in_number = []
        self._barcode_number = 0
       
    def add(self, barcode_label, barcode_sequence):
        if barcode_label in self._barcode_labels:
            raise DuplicateError("Duplicate barcode label '{}'!".
                              format(barcode_label))
        if barcode_sequence in self._barcode_sequences:
            raise DuplicateError("Duplicate barcode sequence '{}'!".
                              format(barcode_sequence))
        if len(barcode_sequence) != self._barcode_length:
            raise ValueError("Length of Barcode '{0}' is not consistent with"\
                             " barcode setting {1}".format(barcode_sequence,
                                                          self._barcode_length))

        self._barcode_sequences.append(barcode_sequence)
        self._barcode_labels.append(barcode_label)
        self._barcode_in_number.append(pattern_to_number(barcode_sequence))
        self._reverse_barcode_in_number.append(
            pattern_to_number(reverse_complement(barcode_sequence)))
        self._barcode_number += 1

    def remove(self, query):
        index_to_be_removed = None
        try:
            index_to_be_removed = self._barcode_labels.index(query)
        except ValueError:
            pass
        try:
            index_to_be_removed = self._barcode_sequences.index(query)
        except ValueError:
            pass
        if index_to_be_removed:
            del self._barcode_labels[index_to_be_removed]
            del self._barcode_sequences[index_to_be_removed]
            del self._barcode_in_number[index_to_be_removed]
            del self._reverse_barcode_in_number[index_to_be_removed]
            return True
        else:
            print("barcode '{}' can not be found in list".format(query))
            return

    def barcode_labels(self):
        return self._barcode_labels

    def barcode_sequences(self):
        return self._barcode_sequences

    def barcode_in_number(self):
        return self._barcode_in_number

    def reverse_barcode_in_number(self):
        return self._reverse_barcode_in_number

    def start_position(self):
        return self._start_position

    def end_position(self):
        return self._end_position

    def reverse_start_position(self):
        return self._reverse_start_position

    def reverse_end_position(self):
        return self._reverse_end_position

    def barcode_number(self):
        return self._barcode_number

    def barcode_set_name(self):
        return self._barcode_set_name

    def insert_size(self):
        return self._insert_size

    def barcode_type(self):
        return self._barcode_type

if __name__ == '__main__':
    main()
