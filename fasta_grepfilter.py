#!/usr/bin/env python
from __future__ import print_function, division
import re
import sys
__author__ = 'markus'


def parse_fasta(fasta):
    while True:
        line = fasta.readline()
        if line[0] == "":
            return
        if line[0] == ">":
            break

    while True:
        header = line[1:].strip("\n")
        seq = ""
        line = fasta.readline()

        while True:
            if not line:
                break
            if line[0] == ">":
                break

            seq += line
            line = fasta.readline()

        yield header, seq.strip("\n")

        if not line:
            return


# and_or(func_list, flag)
#
# Returns a function which goes through
# the functions in func_list in either of
# the following formats:
#    and: (fun1(x) and fun2(x) and fun3(x) ... ) -> bool
#    or:  (fun1(x) or fun2(x) or fun3(x) ... ) -> bool
def and_or(bool_func_list, flag): 
    if flag == 'and':
        bol = False
    elif flag == 'or':
        bol = True
    else:
        raise Exception('Bad flag-argument to function "and_or", must be "and" or "or"')

    def and_or_funcs(x):
        for fun in bool_func_list:
            if fun(x) == bol:
                return bol
        return not bol

    return and_or_funcs

def invert(func): 
    return lambda x: not func(x)


############################
#         FILTERS          #
############################

#  All functions takes the argument(s)
#  for which to use as a threshold
#  and then returns a function which
#  applies that same threshold on its
#  passed Argument and return boolean value.

#
#  Sequence Filters
#

def check_length(min_length): 
    return lambda seq: len(seq) > min_length

def check_substr(substr): 
    substring = re.compile(substr, re.I)
    return lambda seq: bool(substring.search(seq))

def check_content(nucs, min_freq): 
    mf = float(min_freq) / 100

    def content_finder(seq):
        length = len(seq)
        return (len(re.findall(r'[%s]' % nucs, seq, re.I)) / length) > mf
    return content_finder

def check_repeats(repeat_type, repeats):

    # Both args is passed as strings,
    # however, if repeat_type is a string
    # containing only a number, it should be 
    # treated as an int. Consequently it will be
    # interpreted as the definition of the
    # _length_ of the word, not the word itself.
    try:
        repeat_type = int(repeat_type)
    except ValueError:
        pass

    if isinstance(repeat_type, str):
        search_string = '(?P<hit>{0})(?P=hit){1}'.format(repeat_type, '{'+str(repeats)+'}')
    elif isinstance(repeat_type, int):
        search_string = '(?P<hit>[A-Za-z]{0})(?P=hit){1}'.format('{'+str(repeat_type)+'}',
                                                                 '{'+str(repeats)+'}')
    else:
        raise Exception('Bad argument to function "check_repeats", must be string or int')

    rep_finder = re.compile(search_string, re.I)
    return lambda seq: bool(rep_finder.search(seq))

#
# Header-filter
#

def header_filter(pattern_list, mode='', invert_flag=False):

    assert(isinstance(pattern_list, list))

    if mode == 'exact':
        patterns = [re.compile('\b%s\b' % re.escape(x)) for x in pattern_list]
    elif mode == 'regex':
        patterns = [re.compile(x) for x in pattern_list]
    else:
        patterns = [re.compile(re.escape(x)) for x in pattern_list]

    funcs = [lambda x, f=p.search: bool(f(x)) for p in patterns]
    header_pattern_grab = and_or(funcs, 'or')

    if invert_flag: 
        header_pattern_grab = invert(header_pattern_grab)

    return header_pattern_grab

def print_fasta(stats=False):
    # returns a function wich takes two
    # args, the header and the sequence.
    # if stats=True, the returned function
    # will print primitive stats of the sequence,
    # otherwise a normal fasta entry.
    if stats:
        def fasta_printer(head, sequence):
            
            nuc_counts = {'A' : 0,
                          'T' : 0,
                          'G' : 0,
                          'C' : 0}
            ambig = 0
            
            for n in sequence.upper():
                try:
                    nuc_counts[n] += 1
                except KeyError:
                    ambig += 1
            
            print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(head, nuc_counts['A'], nuc_counts['T'],
                                                        nuc_counts['G'], nuc_counts['C'], ambig))
    else:
        def fasta_printer(head, sequence):
            print('>{0}\n{1}'.format(head, sequence))
            
    return fasta_printer

def combine_filters(headerfilter=None, sequencefilter=None):
    
    defaultfunc = lambda x: True
    if headerfilter:
        hf = headerfilter
    else:
        hf = defaultfunc
        
    if sequencefilter:
        sf = sequencefilter       
    else:
        sf = defaultfunc
        
    def combined(header, sequence):
        if hf(header) and sf(sequence):
            return True
        else:
            return False
        
    return combined

if __name__ == '__main__':
    import argparse
    import os.path

    ap = argparse.ArgumentParser()
    ap.add_argument("-f", "--fasta", nargs='?')
    ap.add_argument("-hp", "--header_pattern", metavar='pattern', help="Search for pattern(s) in the header-lines of the file",
                    nargs='+', default='')
    ap.add_argument("-e", "--exact", help="Only match the exact word of HEADER_PATTERN", action="store_true")
    ap.add_argument("-r", "--regex", help="Treat HEADER_PATTERN as a regular expression", action="store_true")
    ap.add_argument("-hf", "--header_file", metavar='rows-of-patterns', help="Pass a list with names to search for in fasta-header.\n"
                                         "Can be combined with -e/-r")
    ap.add_argument("-ss", "--sub_sequence", help="Search for SUB_SEQUENCE in sequences, full regex-support")
    ap.add_argument("-sl", "--sequence_length", help="Search for sequences with a minimum of SEQUENCE_LENGTH", type=int,
                    default=0)
    ap.add_argument("-sc", "--sequence_content", nargs=2, metavar=('nucs', 'percent'), type=str, help="Specify what bases/aminoacids to search for"
                                                                          " and what frequency. I.e -sc gc 50 would return"
                                                                          " entries with a gc-content of at least 50 %%",
                    default='')
                                                                         # Need to write '50 %%' instead
                                                                         # of '50 %' or argparse module will crash.
    ap.add_argument("-sr", "--sequence_repeats", nargs=2, metavar=('word/word-length', 'repetitions'), type=str,
                                                                         help="Search for repetitions of defined word"
                                                                         " or defined length. I. e -sr 3 10 would"
                                                                         " grab entries with 10 trinucleotide repetitions"
                                                                         " and -sr AGG 10 would search for entries with"
                                                                         " AGG repeted 10 times. Note that 'AGGAGGAGG' is TWO"
                                                                         " repetitions of AGG, not three")
    ap.add_argument("-S", "--stats", help="Print tab-delimited statistics", action="store_true")
    ap.add_argument("-I", "--invert", help="Invert the search, so the matched patterns will be discarded instead"
                                           " of those wich isn't matched", action="store_true")
    des = """
    Parse fasta-files and filter them on header-patterns or sequence patterns. Header can be any
    word-pattern, exact word or regex-pattern. Sequence filter can be done on sequence length,
    sub-sequence (i.e. search for specific letter-sequence in each fasta entry), content-frequency
    (ex. gc-content) or repetitions. If option -f/--fasta and a fasta file is omitted, stdin will
    be used instead. Filters can be inverted.
    """
    ap.description = des

    args = ap.parse_args()
    # Locate and retrieve indata
    if args.fasta:
        if not os.path.isfile(args.fasta):
            print('The file %s does not exist, exiting now' % args.fasta, file=sys.stderr)
            exit(1)
        else:
            fasta = open(args.fasta, 'r')
    else:
        fasta = sys.stdin

    # Set header-filter:
    if args.header_file and args.header_pattern:
        print('Cannot handle both a file containing search-patterns and pattern(s) from CLI, choose one.', file=sys.stderr)
        exit(1)

    if args.header_file:
        if not os.path.isfile(args.header_file):
            print('The file %s does not exist, exiting now' % args.file, file=sys.stderr)
            exit(1)

        with open(args.header_file, 'r') as f:
            header_patterns = [x.strip('\n') for x in f.readlines()]

    elif args.header_pattern:
        header_patterns = args.header_pattern
    else:
        header_patterns = None

    if args.exact:
        mode = 'exact'
    elif args.regex:
        mode = 'regex'
    else:
        mode = None

    if header_patterns:
        headerFilter = header_filter(header_patterns, mode, args.invert)
    else:
        headerFilter = None

    # Set sequence-filter:
    seq_filters = []
    if args.sequence_length:
        seq_filters.append(check_length(args.sequence_length))

    if args.sub_sequence:
        seq_filters.append(check_substr(args.sub_sequence))

    if args.sequence_content:
        nucs = args.sequence_content[0]
        freq = args.sequence_content[1]
        seq_filters.append(check_content(nucs, freq))

    if args.sequence_repeats:
        rep_type = args.sequence_repeats[0]
        reps     = args.sequence_repeats[1]
        seq_filters.append(check_repeats(rep_type, reps))

    # Combine sequence filters:
    if seq_filters:
        seqFilter = and_or(seq_filters, 'and')
        if args.invert:
            seqFilter = invert(seqFilter)
    else:
        seqFilter = None

    printer = print_fasta(stats = args.stats)
    combined_filters = combine_filters(headerFilter, seqFilter)
    # Run program:
    if args.stats:
        print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format('name', 'A', 'T', 'G', 'C', 'ambig'))
    for head, seq in parse_fasta(fasta):
        if combined_filters(head, seq):
            printer(head,seq)
    
