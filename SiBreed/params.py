import argparse
import os
import sys
from SiBreed.__init__ import __version__

class Params(object):

    def __init__(self, program_name):
        self.program_name = program_name

    def set_options(self):
        if self.program_name == 'parents':
            parser = self.parents_options()
        elif self.program_name == 'cross':
            parser = self.cross_options()
        elif self.program_name == 'selfing':
            parser = self.selfing_options()
        elif self.program_name == 'backcross':
            parser = self.backcross_options()
        elif self.program_name == 'randcross':
            parser = self.randcross_options()
        elif self.program_name == 'genostat':
            parser = self.genostat_options()
        elif self.program_name == 'genovisual':
            parser = self.genovisual_options()
        elif self.program_name == 'genomap':
            parser = self.genomap_options()
        elif self.program_name == 'phenotype':
            parser = self.phenotype_options()


        if len(sys.argv) == 1:
            args = parser.parse_args(['-h'])
        else:
            args = parser.parse_args()
        return args

    def parents_options(self):
        parser = argparse.ArgumentParser(description='SiBreed version {}'.format(__version__),
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = ('parents -m <TSV file showing linkage map>\n'
                        '        -o <Name of output directory>\n'
                        '        ... \n')

        # set options
        parser.add_argument('-m', '--map',
                            action='store',
                            required=True,
                            type=str,
                            help=('TSV file showing linkage map.\n'
                                  'This file must formatted properly.\n'),
                            metavar='')
        
        parser.add_argument('-o', '--out',
                            action='store',
                            default='result',
                            type=str,
                            help='Directory name of new haploid files.\n',
                            metavar='')
        
        parser.add_argument('-n', '--num',
                            action='store',
                            default=10,
                            type=int,
                            help=('Number of parents generated.\n'
                                  'Default: 10, Max: 35\n'),
                            metavar='')
        
        parser.add_argument('--base_per_character',
                            action='store',
                            default=1000,
                            type=int,
                            help=('Base per one character in haploid files.\n'
                                  'Default: 1000\n'),
                            metavar='')
        
        parser.add_argument('-v', '--version',
                            action='version',
                            version='%(prog)s {}'.format(__version__))
        return parser

    def cross_options(self):
        parser = argparse.ArgumentParser(description='SiBreed version {}'.format(__version__),
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = ('cross -m <TSV file showing linkage map>\n'
                        '      -p1 <File Stem of Parent 1>\n'
                        '      -p2 <File Stem of Parent 2>\n'
                        '      -o <Name of output directory>\n'
                        '      ... \n')

        # set options
        parser.add_argument('-m', '--map',
                            action='store',
                            required=True,
                            type=str,
                            help=('TSV file showing linkage map.\n'
                                  'This file must formatted properly.\n'),
                            metavar='')
        
        parser.add_argument('-p1', '--parent_1',
                            action='store',
                            required=True,
                            type=str,
                            help=('File stem of haploid file of Parent 1.\n'
                                  'For example, if you want to cross haploid file sets named\n'
                                  '"parents/P1_1.hap", "parents/P1_2.hap" and\n'
                                  '"parents/P2_1.hap", "parents/P2_2.hap",\n'
                                  'specify "parents/P1" here.\n'),
                            metavar='')

        parser.add_argument('-p2', '--parent_2',
                            action='store',
                            required=True,
                            type=str,
                            help=('File stem of haploid file of Parent 2.\n'),
                            metavar='')
        
        parser.add_argument('-o', '--out',
                            action='store',
                            default='result',
                            type=str,
                            help='Directory name of new haploid files.\n',
                            metavar='')
        
        parser.add_argument('-s', '--seed',
                            action='store',
                            default=1,
                            type=int,
                            help='Seed number for reproducibility of random numbers.\n',
                            metavar='')
        
        parser.add_argument('-n', '--num',
                            action='store',
                            default=10,
                            type=int,
                            help='Number of children. Default: 10. Max: 999.\n',
                            metavar='')
        
        parser.add_argument('--cpu',
                            action='store',
                            default=1,
                            type=int,
                            help='Number of threads to use\n',
                            metavar='')
        
        
        # set version
        parser.add_argument('-v', '--version',
                            action='version',
                            version='%(prog)s {}'.format(__version__))
        return parser
    
    def selfing_options(self):
        parser = argparse.ArgumentParser(description='SiBreed version {}'.format(__version__),
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = ('selfing -m <TSV file showing linkage map>\n'
                        '        -p <Directory name of haploid file of parents.>\n'
                        '        -o <Name of output directory>\n'
                        '        -n <Number of children>\n'
                        '        ... \n')

        # set options
        parser.add_argument('-m', '--map',
                            action='store',
                            required=True,
                            type=str,
                            help=('TSV file showing linkage map.\n'
                                  'This file must formatted properly.\n'),
                            metavar='')
        
        parser.add_argument('-p', '--parent_directory',
                            action='store',
                            required=True,
                            type=str,
                            help=('Directory name of haploid file of parents.\n'),
                            metavar='')
        
        parser.add_argument('-o', '--out',
                            action='store',
                            default='result',
                            type=str,
                            help='Directory name of new haploid files.\n',
                            metavar='')
        
        parser.add_argument('-s', '--seed',
                            action='store',
                            default=1,
                            type=int,
                            help='Seed number for reproducibility of random numbers.\n',
                            metavar='')
        
        parser.add_argument('-n', '--num',
                            action='store',
                            default=1,
                            type=int,
                            help=('Number of children. Default: 1. Max: 9999.\n'
                                  'Please be careful not to be too much\n'
                                  'when designated directory have many haploid files.\n'),
                            metavar='')
        
        parser.add_argument('--cpu',
                            action='store',
                            default=1,
                            type=int,
                            help='Number of threads to use\n',
                            metavar='')
        
        
        # set version
        parser.add_argument('-v', '--version',
                            action='version',
                            version='%(prog)s {}'.format(__version__))
        return parser

    def backcross_options(self):
        parser = argparse.ArgumentParser(description='SiBreed version {}'.format(__version__),
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = ('backcross -m <TSV file showing linkage map>\n'
                        '          -p <Directory name of donor parents>\n'
                        '          -r <File stem of recurrent parent>\n'
                        '          -o <Name of output directory>\n'
                        '          -n <Number of children>\n'
                        '          ... \n')

        # set options
        parser.add_argument('-m', '--map',
                            action='store',
                            required=True,
                            type=str,
                            help=('TSV file showing linkage map.\n'
                                  'This file must formatted properly.\n'),
                            metavar='')
        
        parser.add_argument('-p', '--parent_directory',
                            action='store',
                            required=True,
                            type=str,
                            help=('Directory name of donor parents.\n'),
                            metavar='')

        parser.add_argument('-r', '--recurrent_parent',
                            action='store',
                            required=True,
                            type=str,
                            help=('File stem of recurrent parent.\n'),
                            metavar='')

        parser.add_argument('-o', '--out',
                            action='store',
                            default='result',
                            type=str,
                            help='Directory name of new haploid files.\n',
                            metavar='')
        
        parser.add_argument('-s', '--seed',
                            action='store',
                            default=1,
                            type=int,
                            help='Seed number for reproducibility of random numbers.\n',
                            metavar='')
        
        parser.add_argument('-n', '--num',
                            action='store',
                            default=1,
                            type=int,
                            help=('Number of children. Default: 1. Max: 999.\n'
                                  'Please be careful not to be too much\n'
                                  'when designated directory have many haploid files.\n'),
                            metavar='')
        
        parser.add_argument('--cpu',
                            action='store',
                            default=1,
                            type=int,
                            help='Number of threads to use\n',
                            metavar='')
        
        
        # set version
        parser.add_argument('-v', '--version',
                            action='version',
                            version='%(prog)s {}'.format(__version__))
        return parser

    def randcross_options(self):
        parser = argparse.ArgumentParser(description='SiBreed version {}'.format(__version__),
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = ('randcross -m <TSV file showing linkage map>\n'
                        '          -pf <File stem of female parent>\n'
                        '          -pm <Directory name of random male parents>\n'
                        '          -o <Name of output directory>\n'
                        '          -n <Number of children>\n'
                        '          ... \n')

        # set options
        parser.add_argument('-m', '--map',
                            action='store',
                            required=True,
                            type=str,
                            help=('TSV file showing linkage map.\n'
                                  'This file must formatted properly.\n'),
                            metavar='')
        
        parser.add_argument('-pf', '--female_parents_directory',
                            action='store',
                            required=True,
                            type=str,
                            help=('Directory name of female parents.\n'),
                            metavar='')

        parser.add_argument('-pm', '--male_parents_directory',
                            action='store',
                            required=True,
                            type=str,
                            help=('Directory name of random male parents.\n'),
                            metavar='')

        parser.add_argument('-o', '--out',
                            action='store',
                            default='result',
                            type=str,
                            help='Directory name of new haploid files.\n',
                            metavar='')
        
        parser.add_argument('-s', '--seed',
                            action='store',
                            default=1,
                            type=int,
                            help='Seed number for reproducibility of random numbers.\n',
                            metavar='')
        
        parser.add_argument('-n', '--num',
                            action='store',
                            default=1,
                            type=int,
                            help=('Number of children. Default: 1. Max: 999.\n'
                                  'Please be careful not to be too much\n'
                                  'when designated directory have many haploid files.\n'),
                            metavar='')
        
        parser.add_argument('--cpu',
                            action='store',
                            default=1,
                            type=int,
                            help='Number of threads to use\n',
                            metavar='')
        
        
        # set version
        parser.add_argument('-v', '--version',
                            action='version',
                            version='%(prog)s {}'.format(__version__))
        return parser

    def genostat_options(self):
        parser = argparse.ArgumentParser(description='SiBreed version {}'.format(__version__),
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = ('genostat -d <Directory containing haploid files>\n'
                        '         ... \n')

        # set options
        parser.add_argument('-d', '--directory',
                            action='store',
                            required=True,
                            type=str,
                            help=('Directory containing haploid files.\n'),
                            metavar='')

        parser.add_argument('--character',
                            action='store',
                            default='0,1',
                            type=str,
                            help=('Characters of haploid files contain.\n'
                                  'Comma seperated value like "0,1" are valid.\n'),
                            metavar='')
        
        # set version
        parser.add_argument('-v', '--version',
                            action='version',
                            version='%(prog)s {}'.format(__version__))
        return parser

    def genovisual_options(self):
        parser = argparse.ArgumentParser(description='SiBreed version {}'.format(__version__),
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = ('genovisual -d <Directory containing haploid files>\n'
                        '           ... \n')

        # set options
        parser.add_argument('-d', '--directory',
                            action='store',
                            required=True,
                            type=str,
                            help=('Directory containing haploid files.\n'),
                            metavar='')
        
        parser.add_argument('--character',
                            action='store',
                            default='0,1',
                            type=str,
                            help=('Characters of haploid files contain.\n'
                                  'Comma seperated value like "0,1" are valid.\n'),
                            metavar='')

        parser.add_argument('--cpu',
                            action='store',
                            default=1,
                            type=int,
                            help='Number of threads to use\n',
                            metavar='')
        
        # set version
        parser.add_argument('-v', '--version',
                            action='version',
                            version='%(prog)s {}'.format(__version__))
        return parser

    def genomap_options(self):
        parser = argparse.ArgumentParser(description='SiBreed version {}'.format(__version__),
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = ('genomap -d <Directory containing haploid files>\n'
                        '         ... \n')

        # set options
        parser.add_argument('-d', '--directory',
                            action='store',
                            required=True,
                            type=str,
                            help=('Directory containing haploid files.\n'),
                            metavar='')

        parser.add_argument('--character',
                            action='store',
                            default='0,1',
                            type=str,
                            help=('Characters of haploid files contain.\n'
                                  'Comma seperated value like "0,1" are valid.\n'),
                            metavar='')
        
        parser.add_argument('-n', '--num',
                            action='store',
                            default=10,
                            type=int,
                            help='Interval of markers (character).\n',
                            metavar='')
        
        # set version
        parser.add_argument('-v', '--version',
                            action='version',
                            version='%(prog)s {}'.format(__version__))
        return parser

    def phenotype_options(self):
        parser = argparse.ArgumentParser(description='SiBreed version {}'.format(__version__),
                                         formatter_class=argparse.RawTextHelpFormatter)
        parser.usage = ('phenotype -q <TSV file showing QTLs>\n'
                        '          -d <Directory containing haploid files>\n'
                        '          ... \n')

        # set options
        parser.add_argument('-q', '--qtl',
                            action='store',
                            required=True,
                            type=str,
                            help=('TSV file showing QTLs.\n'
                                  'This file must formatted properly.\n'),
                            metavar='')
        
        parser.add_argument('-d', '--directory',
                            action='store',
                            required=True,
                            type=str,
                            help=('Directory containing haploid files.\n'),
                            metavar='')
        
        parser.add_argument('--sub',
                            action='store',
                            default='pheno_1',
                            type=str,
                            help='Subscript to identify different analysis.\n',
                            metavar='')

        parser.add_argument('--center',
                            action='store',
                            default=100.0,
                            type=float,
                            help='Standard phenotypic value. default: 100.\n',
                            metavar='')
        
        parser.add_argument('--error',
                            action='store',
                            default=0.0,
                            type=float,
                            help='Variance of error of phenotypic value. default: 0.0.\n',
                            metavar='')

        parser.add_argument('-s', '--seed',
                            action='store',
                            default=1,
                            type=int,
                            help='Seed number for reproducibility of random errors.\n',
                            metavar='')

        parser.add_argument('--cpu',
                            action='store',
                            default=1,
                            type=int,
                            help='Number of threads to use\n',
                            metavar='')
        
        parser.add_argument('-v', '--version',
                            action='version',
                            version='%(prog)s {}'.format(__version__))
        return parser

    def parents_check_args(self, args):
        self.check_file_existance(args.map)
        self.check_dir_absent(args.out)
        if args.num <= 0 or args.num > 36:
            sys.stderr.write('Error. Number of parents generated is invalid.\n')
            sys.exit(1)

    def cross_check_args(self, args):
        self.check_file_existance(args.map)
        self.check_file_existance('{}_1.hap'.format(args.parent_1))
        self.check_file_existance('{}_2.hap'.format(args.parent_1))
        self.check_file_existance('{}_1.hap'.format(args.parent_2))
        self.check_file_existance('{}_2.hap'.format(args.parent_2))
        self.check_dir_absent(args.out)
        if args.num <= 0 or args.num > 999:
            sys.stderr.write('Error. Number of children is invalid.\n')
            sys.exit(1)

    def selfing_check_args(self, args):
        self.check_file_existance(args.map)
        self.check_dir_existance(args.parent_directory)
        self.check_dir_absent(args.out)
        if args.num <= 0 or args.num > 9999:
            sys.stderr.write('Error. Number of children is invalid.\n')
            sys.exit(1)

    def backcross_check_args(self, args):
        self.check_file_existance(args.map)
        self.check_dir_existance(args.parent_directory)
        self.check_file_existance('{}_1.hap'.format(args.recurrent_parent))
        self.check_file_existance('{}_2.hap'.format(args.recurrent_parent))
        self.check_dir_absent(args.out)
        if args.num <= 0 or args.num > 999:
            sys.stderr.write('Error. Number of children is invalid.\n')
            sys.exit(1)

    def randcross_check_args(self, args):
        self.check_file_existance(args.map)
        self.check_dir_existance(args.female_parents_directory)
        self.check_dir_existance(args.male_parents_directory)
        self.check_dir_absent(args.out)
        if args.num <= 0 or args.num > 999:
            sys.stderr.write('Error. Variance of error is invalid.\n')
            sys.exit(1)

    def genostat_check_args(self, args):
        self.check_dir_existance(args.directory)
        self.check_valid_char(args.character)

    def genovisual_check_args(self, args):
        self.check_dir_existance(args.directory)
        self.check_valid_char(args.character)

    def genomap_check_args(self, args):
        self.check_dir_existance(args.directory)
        self.check_valid_char(args.character)

    def phenotype_check_args(self, args):
        self.check_file_existance(args.qtl)
        self.check_dir_existance(args.directory)
        self.check_dir_absent('{}_{}'.format(args.directory, args.sub))
        if args.error < 0:
            sys.stderr.write('Error. Number of children is invalid.\n')
            sys.exit(1)

    def check_file_existance(self, file):
        if not(os.path.isfile(file)):
            sys.stderr.write('Error. Input file does not exist.\n')
            sys.exit(1)

    def check_dir_existance(self, name):
        if not(os.path.isdir(name)):
            sys.stderr.write('Error. Input directory does not exist.\n')
            sys.exit(1)

    def check_dir_absent(self, name):
        if os.path.isdir(name):
            sys.stderr.write('Error. Name of output directory is already used.\n')
            sys.exit(1)

    def check_valid_char(self, input):
        symbols = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
           'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
           'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',
           'u', 'v', 'w', 'y', 'z']
        chars = input.split(',')
        if not set(chars).issubset(symbols):
            sys.stderr.write('Error. Input characters are invalid.\n')
            sys.exit(1)
