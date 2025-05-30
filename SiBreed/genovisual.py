#!/usr/bin/env python3
from multiprocessing import Pool
import os
import sys
from SiBreed import visualization
from SiBreed.params import Params
from SiBreed.utils import time_stamp

pm = Params('genovisual')
args = pm.set_options()

class GenoVisual(object):
    def __init__(self, args):
        self.args = args
        self.dir = args.directory
        self.cpu = args.cpu

        #Check args
        pm.genovisual_check_args(args)

        self.char_1 = args.character.split(',')

        #Make output directory
        self.out = '{}_visual'.format(self.dir)
        os.mkdir(self.out)

    def run(self):
        files = os.listdir(self.dir)
        files_hap_1 = [s for s in files if '_1.hap' in s]
        files_hap_2 = [s for s in files if '_2.hap' in s]
        parent_list = []
        sample_num = 0
        for f in files_hap_1:
            f_stem = '_'.join(f.split('_')[:-1])
            matched = [s for s in files_hap_2 if f_stem + '_2.hap' in s]
            if len(matched) == 1:
                parent_list.append(['{}/{}'.format(self.dir, f_stem), self.char_1, self.out])
                sample_num = sample_num + 1
        if sample_num < 1:
            print(time_stamp(),
                  'Disgnated directory does not cantain valid haploid files.',
                  flush=True)
            sys.exit(1) 

        #Multi processing
        ##############################################
        with Pool(self.cpu) as p:
            p.map(visualization.visualization, parent_list)
        #############################################


def main():
    print(time_stamp(), 'Genotype visualize started.', flush=True)

    prog = GenoVisual(args)
    prog.run()

    print(time_stamp(), 'Successfully finished.\n', flush=True)

if __name__ == '__main__':
    main()
