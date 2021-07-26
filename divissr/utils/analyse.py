#! /usr/bin/env python

# Future
from __future__ import print_function, division

# Generic/Built-in
import sys, os, json
from collections import Counter, defaultdict
import numpy as np

# Owned
try:
    from .seq_utils import kmer_names, get_cycles, build_cycVariations, rev_comp
except:
    from seq_utils import kmer_names, get_cycles, build_cycVariations, rev_comp


def analyse_fasta(args):
    """Generates info JSON """
    repeats_file = args.output
    input_file = args.input

    html_report = os.path.splitext(repeats_file)[0] + '.html'
    print("\nGenerating HTML report. This may take a while..", end="\n\n")

    inf = float('inf')

    repeat_classes = []
    kmer_classes = defaultdict(list)
    cyclical_variations = dict()
    total_repeat_bases = 0
    total_repeat_freq = 0
    longest_lengths = [['seq', 'start', 'stop', 'repeat_class', 0, '+', 0, 'actualrep']]*100
    most_units = [['seq', 'start', 'stop', 'repeat_class', 0, '+', 0, 'actualrep']]*100
    min_length = inf
    min_units = inf

    plot_data = {}
    defaultInfo = {}
    defaultInfo['info'] = { 'SeqInfo': {}, 'RepInfo': {} }

    with open(repeats_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('#'):
                fields = line[1:].split(': ')
                defaultInfo['info']['SeqInfo'][fields[0]] = fields[1]
            else:
                fields = line.split('\t')
                fields = line.split('\t')
                fields[1] = int(fields[1])
                fields[2] = int(fields[2])
                fields[4] = int(fields[4])
                fields[6] = int(fields[6])

                seq_name = fields[0]
                repeat_start = fields[1]
                repeat_end = fields[2]
                repeat_class = fields[3]
                repeat_length = fields[4]
                repeat_ori = fields[5]
                repeat_units = fields[6]
                repeat_actual = fields[7]

                if repeat_class not in repeat_classes:
                    kmer_classes[kmer_names[len(repeat_class)-1]].append(repeat_class)
                    repeat_classes.append(repeat_class)
                    cyclical_variations[repeat_class] = build_cycVariations(repeat_class)

                    plot_data[repeat_class] = dict()
                    plot_data[repeat_class][repeat_length] = [0]*len(cyclical_variations[repeat_class])

                total_repeat_bases += repeat_length
                total_repeat_freq += 1

                if repeat_length not in plot_data[repeat_class]: 
                    plot_data[repeat_class][repeat_length] = \
                        [0]*len(cyclical_variations[repeat_class])

                plot_data[repeat_class][repeat_length] \
                    [cyclical_variations[repeat_class].index(repeat_actual)] += 1

                if min_units > repeat_units: min_units = repeat_units
                if min_length > repeat_length: min_length = repeat_length

                if (longest_lengths[-1][4] < repeat_length) or \
                (longest_lengths[-1][4] == repeat_length and \
                repeat_class < longest_lengths[-1][3]):
                    longest_lengths[-1] = fields
                    longest_lengths.sort(key=lambda x: x[4])
                    longest_lengths.reverse()
                
                if (most_units[-1][6] < repeat_units) or \
                (most_units[-1][6] == repeat_units and \
                repeat_class < longest_lengths[-1][3]):
                    most_units[-1] = fields
                    most_units.sort(key=lambda x: x[6])
                    most_units.reverse()
    repeat_options = ""
    kmers = list(kmer_classes.keys())
    kmers.sort(key=lambda x: kmer_names.index(x))
    for kmer in kmers:
        repeat_options += '<optgroup label="%s">' %(kmer)
        for r in kmer_classes[kmer]:
            repeat_options += '<option value="%s">%s</option>' %(r, r)
        repeat_options += '</optgroup>'

    total_bases = int(defaultInfo['info']['SeqInfo']['GenomeSize'])
    defaultInfo['info']['RepInfo']['PlotData'] = plot_data
    defaultInfo['info']['RepInfo']['NumRepClasses'] = len(plot_data.keys())
    defaultInfo['info']['RepInfo']['TotalRepBases'] = total_repeat_bases
    defaultInfo['info']['RepInfo']['TotalRepFreq'] = total_repeat_freq
    defaultInfo['info']['RepInfo']['PercentGenomeCovered'] = \
      str(round((total_repeat_bases/total_bases)*100, 2)) + "%"
    defaultInfo['info']['RepInfo']['RepDensityByFreq'] = \
      round((total_repeat_freq/total_bases)*1000000, 2)
    defaultInfo['info']['RepInfo']['RepDensityByBases'] = \
      round((total_repeat_bases/total_bases)*1000000, 2)
    defaultInfo['info']['RepInfo']['MinLength'] = min_length
    defaultInfo['info']['RepInfo']['MinUnits'] = min_units
    defaultInfo['info']['RepInfo']['LongestRepeats'] = []
    defaultInfo['info']['RepInfo']['MostRepeatUnits'] = []
    defaultInfo['info']['RepInfo']['AllRepClasses'] = repeat_classes
    for a in longest_lengths:
        testDict = {
            'seq': a[0], 'start': a[1], 'stop': a[2], 'repClass': a[3], 
            'repLength': a[4], 'repOri': a[5], 'repUnit': a[6], 'actualRep': a[7]
        }
        defaultInfo['info']['RepInfo']['LongestRepeats'].append(testDict)
    for a in most_units:
        testDict = {
            'seq': a[0], 'start': a[1], 'stop': a[2], 'repClass': a[3], 
            'repLength': a[4], 'repOri': a[5], 'repUnit': a[6], 'actualRep': a[7]
        }
        defaultInfo['info']['RepInfo']['MostRepeatUnits'].append(testDict)
    defaultInfo = 'const data =' + json.dumps(defaultInfo)
    writetoHTML(html_report, defaultInfo, repeat_options, 'fasta')


def analyse_fastq(args):

    """Generates HTML report for fastq files."""
    input_file = args.input
    repeats_file = args.output

    defaultInfo = {'info': { 'SeqInfo': {}, 'RepInfo': {} }}

    html_report = os.path.splitext(repeats_file)[0] + '.html'
    print("\nGenerating HTML report. This may take a while..", end="\n\n")

    inf = float('inf')
    min_length = inf
    max_length = 0
    
    repeat_classes = []
    kmer_classes = defaultdict(list)
    defaultInfo['info']['SeqInfo']['FileName'] = args.input.split('/')[-1]
    plot_data = {}
    with open(repeats_file, 'r') as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('#'):
                fields = line[1:].split(': ')
                if 'Rep' in fields[0]:
                    defaultInfo['info']['RepInfo'][fields[0]] = fields[1]
                else:
                    defaultInfo['info']['SeqInfo'][fields[0]] = fields[1]
            else:
                fields = line.strip().split('\t')
                rclass = fields[0]
                rfreq = int(fields[1]) 
                rreads = int(fields[2]) 
                rbases = int(fields[3])
                repeat_classes.append(rclass)
                plot_data[rclass] = { 'Freq': rfreq, 'Reads': rreads, 'Bases': rbases }
                kmer_classes[kmer_names[len(rclass)-1]].append(rclass)
                length_freq = {}
                for l in fields[-1].split(';'):
                    length = int(l.split('-')[0])
                    if length < min_length: min_length = length
                    if length > max_length: max_length = length
                    freq = int(l.split('-')[1])
                    length_freq[length] = freq
                plot_data[rclass]['LenFreq'] = length_freq
    
    defaultInfo['info']['RepInfo']['MinLength'] = min_length
    defaultInfo['info']['RepInfo']['MaxLength'] = max_length
    defaultInfo['info']['RepInfo']['PlotData'] = plot_data
    defaultInfo['info']['RepInfo']['AllRepClasses'] = repeat_classes
    read_length_dist_text = defaultInfo['info']['SeqInfo']['ReadLengthDist']  # removing the last semicolon
    read_length_dist = {}
    for l in read_length_dist_text.split(';'):
        read_length_dist[int(l.split('-')[0])] = int(l.split('-')[1])

    defaultInfo['info']['SeqInfo']['ReadLengthDist'] = read_length_dist
    repeat_options = ""
    kmers = list(kmer_classes.keys())
    kmers.sort(key=lambda x: kmer_names.index(x))
    for kmer in kmers:
        repeat_options += '<optgroup label="%s">' %(kmer)
        for r in kmer_classes[kmer]:
            repeat_options += '<option value="%s">%s</option>' %(r, r)
        repeat_options += '</optgroup>'
    defaultInfo = 'const data =' + json.dumps(defaultInfo)
    writetoHTML(html_report, defaultInfo, repeat_options, 'fastq')


def writetoHTML(html_file, defaultInfo, repeat_options, fformat):
    """Collates the scripts needed for the template html file."""
    
    html_handle = open(html_file, 'w')
    current_dir = os.path.dirname(__file__)
    
    template = open('%s/../lib/looper_%s_template.html' %(current_dir, fformat), 'r').read()
    fontawesome_js = open('%s/../lib/src/all.js' %(current_dir), 'r').read()
    multiselect_css = open('%s/../lib/styles/multi-select.min.css' %(current_dir), 'r').read()
    main_css = open('%s/../lib/styles/main.css' %(current_dir), 'r').read()

    jquery_js = open("%s/../lib/src/jquery-3.5.0.min.js" %(current_dir), "r").read()
    multiselect_js = open('%s/../lib/src/jquery.multi-select.min.js' %(current_dir), 'r').read()
    echarts_js = open('%s/../lib/src/echarts_cp1252.min.js' %(current_dir), 'r').read()
    main_js = open('%s/../lib/src/main_%s.js' %(current_dir, fformat), 'r').read()
    tables_js = open('%s/../lib/src/tables_%s.js' %(current_dir, fformat), 'r').read()

    template = template.format(
        fontawesome_js = fontawesome_js,
        multiselect_css = multiselect_css, 
        main_css = main_css, 
        jquery_js = jquery_js, 
        multiselect_js = multiselect_js, 
        echarts_js = echarts_js, 
        analyse_data_js = defaultInfo, 
        main_js = main_js, 
        tables_js = tables_js, 
        repeat_options = repeat_options,
    )

    print(template, file=html_handle)
    html_handle.close()
    print("HTML report successfully saved to " + html_file)
