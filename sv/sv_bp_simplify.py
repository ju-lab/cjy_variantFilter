'''this is for using svclone input in the 'annotate' step
'''
import os
import re
import argparse

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

    
    @classmethod
    def fromstring(cls, position_string):
        chromosome = position_string.split(':')[0]
        start = int(position_string.split(':')[1].split('-')[0])
        end = int(position_string.split(':')[1].split('-')[1])
        return Position(chromosome, start, end)



def write_simple_bp(inputfile, output_dir):

    outfile_base = re.sub(r'.tsv', '.simplebp.tsv', os.path.basename(inputfile))
    outfile = os.path.join(output_dir, outfile_base)
    print(outfile)
    with open(outfile, 'w') as g:
        g.write('chr1\tpos1\tchr2\tpos2\n')
        with open(inputfile, 'r') as f:
            for line in f:
                print(line);bp1, bp2, svtype, sampleName = line.strip().split()
                bp1Pos = Position.fromstring(bp1)
                bp2Pos = Position.fromstring(bp2)
                g.write(str(bp1Pos.chromosome) + '\t' + str(bp1Pos.start) + '\t' + str(bp2Pos.chromosome) + '\t' + str(bp2Pos.start) + '\n')


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputfile', required=True)
    parser.add_argument('-o', '--outputdir', required=False, default=os.getcwd())

    args = vars(parser.parse_args())

    return args['inputfile'], args['outputdir']

def main():
    inputfile, output_dir = argument_parser()
    print(output_dir) 
    write_simple_bp(inputfile, output_dir)


if __name__=='__main__':
    main()


