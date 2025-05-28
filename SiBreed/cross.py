#!/usr/bin/env python3
from multiprocessing import Pool
import os
from SiBreed import recombination, rwhap
from SiBreed.params import Params
from SiBreed.utils import command, time_stamp

pm = Params('cross')
args = pm.set_options()

class Cross(object):
    def __init__(self, args):
        self.args = args
        self.p1 = args.parent_1
        self.p2 = args.parent_2
        self.out = args.out
        self.seed = args.seed
        self.num = args.num
        self.cpu = args.cpu
        
        #Check args
        pm.cross_check_args(args)

        #Read map file
        self.map = rwhap.readmap(args.map)

        #Make output directory
        os.mkdir('{}'.format(self.out))

        #Write command infomation
        command(self.out)

    def run(self):
        p1_stem = self.p1.split('/')[-1] #dir/P1 --> P1
        p2_stem = self.p2.split('/')[-1]
        throw = []
        for i in range(self.num):
            if self.num < 10:
                num_str = '{:01}'.format(i+1)
            elif self.num < 100:
                num_str = '{:02}'.format(i+1)
            else:
                num_str = '{:03}'.format(i+1)
            stem = '{}x{}-{}'.format(p1_stem, p2_stem, num_str)
            one = [self.map, self.p1, self.p2, stem, self.seed, self.out]
            throw.append(one)

        #Multi processing
        ##############################################
        with Pool(self.cpu) as p:
            p.map(recombination.recombination, throw)
        #############################################

def main():
    print(time_stamp(), 'Cross simulation started.', flush=True)

    prog = Cross(args)
    prog.run()

    print(time_stamp(), 'Successfully finished.\n', flush=True)

if __name__ == '__main__':
    main()
