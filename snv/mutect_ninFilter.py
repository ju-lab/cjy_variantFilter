"""
March 12 2018 Jongsoo Yoon
Mutect Filter, that looks at a presence of variant reads in the normal
"""

import argparse
import cyvcf2
import re
import time
import datetime
import os
import sys
import vcf



def argument_parser():
    """parses argument passed on from command line"""
    parser = argparse.ArgumentParser(
        description='separate mutect only calls, strelka only calls, and intersect calls from two vcf files')

    parser.add_argument('-o', '--outputDIR', required=False, default=os.getcwd())
    parser.add_argument('-m', '--mutect', required=True, help='Mutect VCF')
    parser.add_argument('-n', '--normal_column', type=int, required=True, help='Tells whether normal column is written before cancer column. For example after the <FORMAT> column of vcf, if the order is <TUMOR> <NORMAL>, then --normal_column is 1. If <NORMAL> <TUMOR> then --normal_column is 0')
    parser.add_argument('-t', '--threshold', type=float, required=True, help ='Threshold to filter out those variants that has variant reads in normal Bam')

    args = vars(parser.parse_args())

    outputDIR = args['outputDIR']
    mutect = args['mutect']
    normal_column = args['normal_column']
    threshold = args['threshold']
    return outputDIR, mutect, normal_column, threshold

def filter_normal(mutect_vcf, output_dir,  normal_column, threshold=0.01):
    # nin filter = Not In Normal Filter 
    filtered_mutect = os.path.join(output_dir, re.sub(string=os.path.basename(mutect_vcf), pattern=r'.vcf$', repl='.nin_filter_' + str(threshold) + '.vcf'))
    print(filtered_mutect)
    # copy vcf header from input vcf
    os.system('bcftools view -h ' + mutect_vcf+ ' > ' + filtered_mutect)

    with open(filtered_mutect, 'a') as h:
        for variant in cyvcf2.VCF(mutect_vcf):
            normal_vaf = variant.format('FA')[normal_column]
            if normal_vaf > threshold:
                # normal VAF is too high, filter out
                print('filtered: ' + str(variant)  + str(normal_vaf))
                
            else:
#                print(str(variant)  + str(normal_vaf))
                h.write(str(variant))

    return 0

def main():
    outputDIR, mutect_vcf, normal_column, threshold = argument_parser()
    filter_normal(mutect_vcf, outputDIR,  normal_column, threshold)
    return 0

if __name__=='__main__':
    main()

