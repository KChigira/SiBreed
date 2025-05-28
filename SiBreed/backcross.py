#!/usr/bin/env python3
from multiprocessing import Pool
import os
from SiBreed import recombination, rwhap
from SiBreed.params import Params
from SiBreed.utils import command, time_stamp

pm = Params('backcross')
args = pm.set_options()

class Backcross(object):
    def __init__(self, args):
        self.args = args
        self.dir = args.parent_directory
        self.p2 = args.recurrent_parent
        self.out = args.out
        self.seed = args.seed
        self.num = args.num
        self.cpu = args.cpu
        
        #Check args
        pm.backcross_check_args(args)

        #Read map file
        self.map = rwhap.readmap(args.map)

        #Make output directory
        os.mkdir('{}'.format(self.out))

        #Write command infomation
        command(self.out)

    def run(self):
        files = os.listdir(self.dir)
        files_hap_1 = [s for s in files if '_1.hap' in s]
        files_hap_2 = [s for s in files if '_2.hap' in s]
        parent_list = []
        stem_list = []
        for f in files_hap_1:
            f_stem = ''.join(f.split('_')[:-1])
            matched = [s for s in files_hap_2 if f_stem + '_2.hap' in s]
            if len(matched) == 1:
                parent_list.append('{}/{}'.format(self.dir, f_stem))
                stem_list.append(f_stem)
        p2_stem = self.p2.split('/')[-1]

        throw = []
        for h in range(len(parent_list)):
            for i in range(self.num):
                if self.num < 10:
                    num_str = '{:01}'.format(i+1)
                elif self.num < 100:
                    num_str = '{:02}'.format(i+1)
                else:
                    num_str = '{:03}'.format(i+1)
                stem = '{}x(bc){}-{}'.format(stem_list[h], p2_stem, num_str)
                one = [self.map, parent_list[h], self.p2, 
                       stem, self.seed, self.out]
                throw.append(one)

        #Multi processing
        ##############################################
        with Pool(self.cpu) as p:
            p.map(recombination.recombination, throw)
        #############################################

def main():
    print(time_stamp(), 'Backcross simulation started.', flush=True)

    prog = Backcross(args)
    prog.run()

    print(time_stamp(), 'Successfully finished.\n', flush=True)

if __name__ == '__main__':
    main()
