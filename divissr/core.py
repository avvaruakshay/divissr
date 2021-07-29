#!/usr/bin/env python

"""
@author: Avvaru Akshay
@filename: core.py

DiviSSR python script.

This code calculates and inserts the desired parameters in C++ code as constants.
Thus the C++ code doesn't lag on calculating the parameters during the runtime
and also optimises the divisions as the divisors are known at the compile time.
Compiles the generated c++ code and runs repeat indentification using the
compiled c++ code on input sequence file.

"""

# Future
from __future__ import print_function, division

# Generic/Built-in
import os, argparse, sys
from os.path import splitext
from datetime import datetime

# Owned
try:
    from .utils.analyse import analyse_fasta, analyse_fastq
    from .utils.annotation import annotate_repeats
    from .utils.seq_utils import rawcharCount, getGenomeInfo
except:
    from utils.analyse import analyse_fasta, analyse_fastq
    from utils.annotation import annotate_repeats
    from utils.seq_utils import rawcharCount, getGenomeInfo

if sys.version_info[0] == 2:
    pass


def motif_factors(m, M):
    """
    Generates the set of non redundant motif sizes to be checked for division
    rule.

    Parameters
    ----------
    m : <int> minimum motif size

    M : <int> maximum motif size

    Returns
    -------
    factors : <list> sorted(descending) list of non-redundant motifs.
    """
    factors = []
    for i in range(M, m-1, -1):
        check = 0
        for j in factors:
            if j % i == 0: check = 1; break
        if check == 0: factors.append(i)
    
    return factors

def getArgs():
    """
    Parses command line arguments and returns them to the caller
    """
    __version__ = 'v0.1.1'
    parser = argparse.ArgumentParser()
    parser._action_groups.pop()

    required = parser.add_argument_group('Required arguments')
    required.add_argument('-i', '--input', required=True, metavar='<FILE>', help='Input sequence file.')
    
    optional = parser.add_argument_group('Optional arguments')
    
    #Basic options
    optional.add_argument('-o', '--output', metavar='<FILE>', help='Output file name. Default: Input file name + _divissr.tsv')
    optional.add_argument('--format', metavar='<STR>', default='fasta', help='Input file format. Default: fasta, Permissible: fasta, fastq')
    optional.add_argument('-v-', '--version', action='version', version='divissr ' + __version__)
        
    #Selections options based on motif size and seq lengths
    optional.add_argument('-m', '--min-motif-size', type=int, metavar='<INT>', default=1, help='Minimum size of a repeat motif in bp (Not allowed with -rep)')
    optional.add_argument('-M', '--max-motif-size', type=int, metavar='<INT>', default=6, help='Maximum size of a repeat motif in bp (Not allowed with -rep)')
    optional.add_argument('-l', '--min-length', type=int, metavar='<INT>', help='Minimum length cutoff of repeat')
    optional.add_argument('--filter-reads', action='store_true', default=False, help='Seperate out the reads with repeats in a new file.')

    #coumpound repeat
    compound = parser.add_argument_group('Compound repeat arguments')
    compound.add_argument('--compound', action='store_true', default=False, help='Report compound repeats. The output of compound repeats is in a separate file with the suffix ".compound".')
    compound.add_argument('-d', '--comp-dist', type=int, metavar='<INT>', default=0, help='Maximum distance between individual repeats of compound repeat. Use negative to denote overlap. Default: 0')
    
    # Analysis options
    optional.add_argument('-a', '--analyse', action='store_true', default=False, help='Generate a summary HTML report.')

    # Annotation options
    annotation = parser.add_argument_group('Annotation arguments')
    annotation.add_argument('-g', '--annotate', metavar='<FILE>', help='Genic annotation input file for annotation, Both GFF and GTF can be processed. Use --anno-format to specify format.')
    annotation.add_argument('--anno-format', metavar='<STR>',default='GFF', type=str, help='Format of genic annotation file. Valid inputs: GFF, GTF. Default: GFF')
    annotation.add_argument('--gene-key', metavar='<STR>', default=None, type=str, help='Attribute key for geneId. The default identifier is "gene" for GFF and "gene_id" for GTF. Please check the annotation file and pick a robust gene identifier from the attribute column.')
    annotation.add_argument('--up-promoter', metavar='<INT>', type=int, default=1000, help='Upstream distance(bp) from TSS to be considered as promoter region. Default 1000')
    annotation.add_argument('--down-promoter', metavar='<INT>', type=int, default=1000, help='Downstream distance(bp) from TSS to be considered as promoter region. Default 1000')

    args = parser.parse_args()

    if args.min_length is None:
        args.min_length = 2 * args.max_motif_size
    
    if args.output is None:
        args.output = splitext(args.input)[0] + '_divissr.tsv'

    if args.gene_key is None:
        if args.anno_format == 'GFF': args.gene_key = 'gene'
        elif args.anno_format == 'GTF': args.gene_key = 'gene_id'
    
    print('Printing output to:', args.output)

    return args

def main():
    """Main function of divissr"""

    start_time = datetime.now()
    current_dir = os.path.dirname(__file__)
    current_dir = '/'.join(os.path.abspath(__file__).split('/')[:-1])

    args = getArgs()

    m = args.min_motif_size
    M = args.max_motif_size
    cutoff = args.min_length

    if sys.platform.startswith('win'):
        current_dir = '\\'.join(os.path.abspath(__file__).split('\\')[:-1])
        if args.format == 'fasta':
            command = ['{current_dir}\\templates\\wdivissr_fasta.exe'.format(current_dir=current_dir)]
            if args.analyse: command.append('-a')
            if args.compound:
                command.append('--compound')
                comp_out = '.'.join(args.output.split('.')[:-1]) + '.compound.' + args.output.split('.')[-1]
                command.append('--comp-out {comp_out}'.format(comp_out=comp_out))
                command.append('--comp-dist {comp_dist}'.format(comp_dist=args.comp_dist))
        if args.format == 'fastq':
            command = ['{current_dir}\\templates\\wdivissr_fastq.exe'.format(current_dir=current_dir)]
        command.append('-i {input_file}'.format(input_file=args.input))
        command.append('-o {output_file}'.format(output_file=args.output))
        command.append('-m {m}'.format(m=m))
        command.append('-M {M}'.format(M=M))
        command.append('-l {l}'.format(l=cutoff))
        os.system(' '.join(command))

    
    elif sys.platform.startswith('linux'):
        motif_checks = motif_factors(m, M)
        N = len(motif_checks)
        motif_checks_string = '{'
        for a in motif_checks: motif_checks_string += str(a) + ','
        motif_checks_string = motif_checks_string[:-1] + '}'

        divisor = []
        rem_shift = []
        # Calculating all division variables based on the input parameters
        for i, a in enumerate(motif_checks): 
            d = cutoff // motif_checks[i]
            r = cutoff % motif_checks[i]
            D = ('0'*((2*a)-1) + '1')*d + '0'*(2*r)
            divisor.append(int(D, 2))
            rem_shift.append( int(2*(cutoff - (cutoff % motif_checks[i]))) )

        compound_string = 'uint compound = 0;'
        overlap_d_string = 'int overlap_d = 0;'
        if args.compound:
            compound_string = compound_string[:-2] + '1;'
            overlap_d_string = overlap_d_string[:-2] + str(args.comp_dist) + ';'

        div_string = 'uint64_t divisor[{N}] = '.format(N=N) + '{'
        for a in divisor: div_string += str(a) + ','
        div_string = div_string[:-1] + '}'

        rem_shift_string = 'uint rem_shift[{N}] = '.format(N=N) + '{'
        for a in rem_shift: rem_shift_string += str(a) + ','
        rem_shift_string = rem_shift_string[:-1] + '}'

        template_file = ''
        filtered_file = ''
        seq_count = 0
        if args.format == 'fasta':
            seq_count = rawcharCount(args.input, '>')
            if args.input.endswith('gz'):
                template_file = '{current_dir}/templates/pydivissr_fasta_gzip_template.cpp'.format(current_dir=current_dir)
            else:
                template_file = '{current_dir}/templates/pydivissr_fasta_template.cpp'.format(current_dir=current_dir)
        elif args.format == 'fastq':
            if args.input.endswith('.gz'):
                template_file = '{current_dir}/templates/pydivissr_fastq_gzip_template.cpp'.format(current_dir=current_dir)
            else:
                template_file = '{current_dir}/templates/pydivissr_fastq_template.cpp'.format(current_dir=current_dir)
            if args.filter_reads: 
                filtered_file = splitext(args.input)[0] + '_divissr.filtered.fastq'

        script = open('{current_dir}/pydivissr.cpp'.format(current_dir=current_dir), 'w')
        with open(template_file) as fh:
            for line in fh:
                line = line.rstrip()
                if line.strip() == '$ python_input;':
                    prefix = line[:line.find('$')]
                    line = '{prefix}uint cutoff = {cutoff};'.format(prefix=prefix, cutoff=cutoff)
                    line += '\n{prefix}uint m = {m};'.format(prefix=prefix, m=m)
                    line += '\n{prefix}uint M = {M};'.format(prefix=prefix, M=M)
                    line += '\n{prefix}uint motif_checks[{N}] = {motif_checks_string};'.format(prefix=prefix, N=N, motif_checks_string=motif_checks_string)
                    line += '\n{prefix}uint N = {N};'.format(prefix=prefix, N=N)
                    line += '\n{prefix}{div_string};'.format(prefix=prefix, div_string=div_string)
                    line += '\n{prefix}{rem_shift_string};'.format(prefix=prefix, rem_shift_string=rem_shift_string)
                    line += '\n{prefix}{compound_string};'.format(prefix=prefix, compound_string=compound_string)
                    line += '\n{prefix}{overlap_d_string};'.format(prefix=prefix, overlap_d_string=overlap_d_string)
                    line += '\n{prefix}int sequences = {seq_count};'.format(prefix=prefix, seq_count=seq_count)
                if line.strip() == '$ fasta_gzip;':
                    prefix = line[:line.find('$')]
                    line = '{prefix}string fin = "{filename}";'.format(prefix=prefix, filename=args.input)
                print(line, file=script)
        script.close()
        status = os.system('g++ {current_dir}/pydivissr.cpp -O3 -o {current_dir}/pydivissr'.format(current_dir=current_dir))
        if sys.platform.startswith('linux') and status != 0:
            print('\n\nError compiling the cpp auxiliary code. Please check for proper installation of g++. ')
            sys.exit()
        
        end_time = datetime.now()
        total_time = end_time-start_time

        print('\nProcessing time: {total_time} secs'.format(total_time=round(total_time.total_seconds(), 2)))

        if args.format == 'fasta':
            analyse_flag = 0
            if args.analyse:analyse_flag = 1
            comp_out = ''
            if args.compound:
                comp_out = '.'.join(args.output.split('.')[:-1]) + '.compound.' + args.output.split('.')[-1]
            if args.input.endswith('.gz'):
                os.system(
                    'zcat {input} | {current_dir}/pydivissr {output} {analyse_flag} {comp_out}'
                        .format(
                            input=args.input,
                            current_dir=current_dir,
                            output=args.output,
                            analyse_flag=analyse_flag,
                            comp_out=comp_out
                        )
                )
            else:
                os.system(
                    '{current_dir}/pydivissr {input} {output} {analyse_flag} {comp_out}'
                        .format(
                            current_dir=current_dir,
                            input=args.input,
                            output=args.output,
                            analyse_flag=analyse_flag,
                            comp_out=comp_out
                        )
                )

        elif args.format == 'fastq':
            if args.input.endswith('.gz') and args.filter_reads:
                os.system('zcat {input} {filtered_file} | {current_dir}/pydivissr {output} '.format(input=args.input, current_dir=current_dir, output=args.output, filtered_file=filtered_file))
            elif args.input.endswith('.gz'):
                os.system('zcat {input} | {current_dir}/pydivissr {output} '.format(input=args.input, current_dir=current_dir, output=args.output))
            elif args.filter_reads:
                os.system('{current_dir}/pydivissr  {input} {output} {filtered_file}'.format(current_dir=current_dir, input=args.input, output=args.output, filtered_file=filtered_file))
            else:
                os.system('{current_dir}/pydivissr {input} {output}'.format(current_dir=current_dir, input=args.input, output=args.output))

    if args.format == 'fasta' and args.annotate: annotate_repeats(args)

    if args.analyse:
        if args.format == 'fasta': analyse_fasta(args)
        elif args.format == 'fastq': analyse_fastq(args)
    
    end_time = datetime.now()
    total_time = end_time-start_time    
    print('\nTotal time taken: {total_time}'.format(total_time=total_time))


if __name__ == "__main__":
    main()
