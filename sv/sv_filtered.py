'''This script writes the filtered SV output in the original VCF format of a given file. 
Checks if the SV is present in the filtered_tsv which is a tabulated output of the intersect between the
original VCF input and another VCF from a different SV caller. 
Output only consists of those SVs that are present from both callers
'''
import cyvcf2
import argparse
import re
import os


class Position():
    ''' python class for handling genomic positions
    0-based
    '''
    def __init__(self, chromosome, start, end, is_bp=None, clipped_reads=None):
        self.chromosome = chromosome
        self.start = start
        self.end = end

    def __repr__(self):
        return str(self.chromosome) + ":" + str(self.start) + '-' + str(self.end)

    def __str__(self):
        return str(self.chromosome) + ":" + str(self.start) + '-' + str(self.end)
   
    def __hash__(self):
        return hash((self.chromosome, self.start, self.end))

    def __eq__(self, other):
        if isinstance(other, Position):
            return self.chromosome == other.chromosome and self.start == other.start and self.end == other.end
        else:
            print("Not of the same class, cannot compare equality")
            return None

    def __lt__(self, other):
        if isinstance(other, Position):
            if self.chromosome < other.chromosome:
                return True
            elif self.chromosome == other.chromosome:
                if self.start < other.start:
                    return True
                else:
                    return False
            else:
                return False

    
    @classmethod
    def fromstring(cls, position_string):
        chromosome = position_string.split(':')[0]
        start = int(position_string.split(':')[1].split('-')[0])
        end = int(position_string.split(':')[1].split('-')[1])
        return Position(chromosome, start, end)

    @staticmethod
    def overlap(position1, position2):
        '''true if position1 and position2 has more than 1 overlapping base'''
        try:
            if isinstance(position1, Position) and isinstance(position2, Position):
                if position1.chromosome == position2.chromosome:
                    if min(position1.end, position2.end) > max (position1.start, position2.start):
                        return True
                    else:
                        return False
                else:
                    return False  # cannot compare if two positions are in different chromosome
            else:
                return None # has to be Posiiton class.
        except:
            Exception
    
    def extend(self, direction, basepairs):
        """extends objects in by specified base pairs, either upstream, downstream, or both"""
        if direction=="up":
            return Position(self.chromosome, max(0, self.start-basepairs), end)
        elif direction=="down":
            return Position(self.chromosome, self.start, self.end + basepairs)
        elif direction=="both":
            return Position(self.chromosome, max(0, self.start - basepairs), self.end + basepairs)
        else:
            print('direction has to be either up, down, or both')
            raise ValueError



def argument_parser():
    """parses argument passed on from command line"""
    parser = argparse.ArgumentParser(
        description='Find fusion evidence from WGS results')

    parser.add_argument('-v', '--vcf', required=True, help='VCF file to filter')
    parser.add_argument('-o', '--outputDIR', required=False, default=os.getcwd())

    parser.add_argument('-f', '--filtered_tsv', required=True, help='Filtered tsv file that has [breakpoint1]\t[breakpoint2]\t[svtype]\t[sampleName] format')
    args = vars(parser.parse_args())

    vcf = args['vcf']
    outputDIR = args['outputDIR']
    filtered_tsv = args['filtered_tsv']
    return vcf, outputDIR, filtered_tsv


def filtered_breakpoints(filtered_tsv):
    bp1_set = set()
    with open(filtered_tsv, 'r') as f:
        for line in f:
            bp1, bp2, svtype, sampleName = line.strip().split()
            bp1Pos = Position.fromstring(bp1)
            bp2Pos = Position.fromstring(bp2)

            bp1_set.add(bp1Pos)
    
    return bp1_set


def write_filtered_vcf(original_vcf, filtered_tsv, outputDIR):
    filtered_bp = filtered_breakpoints(filtered_tsv)
    filtered_output_basename = re.sub(string=os.path.basename(original_vcf), pattern=r'.vcf$|.vcf.gz$', repl='.intersectFiltered.vcf')
    filtered_output = os.path.join(outputDIR, filtered_output_basename)
    
    # create writer object for the output
    filtered_vcf = cyvcf2.Writer(filtered_output, cyvcf2.VCF(original_vcf))

    #
    for variant in cyvcf2.VCF(original_vcf):
        variant_bp1 = Position(variant.CHROM, variant.POS, variant.POS + 1)
        if variant_bp1 in filtered_bp:
            filtered_vcf.write_record(variant)
        else:
            pass
            # these variants do not have support from the 2nd SV vcf file, thus discarded from the output

def main():
    original_vcf, outputDIR, filtered_tsv = argument_parser()
    write_filtered_vcf(original_vcf, filtered_tsv, outputDIR)

if __name__ == '__main__':
    main()
    


