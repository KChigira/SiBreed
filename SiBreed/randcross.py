#!/usr/bin/env python3
import hashlib
from multiprocessing import Pool
import os
import random
from SiBreed import recombination, rwhap
from SiBreed.params import Params
from SiBreed.utils import command, time_stamp

pm = Params('randcross')
args = pm.set_options()

class Randcross(object):
    def __init__(self, args):
        self.args = args
        self.fem_dir = args.female_parents_directory
        self.mal_dir = args.male_parents_directory
        self.out = args.out
        self.seed = args.seed
        self.num = args.num
        self.cpu = args.cpu
        
        #Check args
        pm.randcross_check_args(args)

        #Read map file
        self.map = rwhap.readmap(args.map)

        #Make output directory
        os.mkdir('{}'.format(self.out))

        #Write command infomation
        command(self.out)

    def run(self):
        #Female parents
        files = os.listdir(self.fem_dir)
        files_hap_1 = [s for s in files if '_1.hap' in s]
        files_hap_2 = [s for s in files if '_2.hap' in s]
        fem_parent_list = []
        fem_stem_list = []
        for f in files_hap_1:
            f_stem = ''.join(f.split('_')[:-1])
            matched = [s for s in files_hap_2 if f_stem + '_2.hap' in s]
            if len(matched) == 1:
                fem_parent_list.append('{}/{}'.format(self.fem_dir, f_stem))
                fem_stem_list.append(f_stem)
        
        #Male parents
        files = os.listdir(self.mal_dir)
        files_hap_1 = [s for s in files if '_1.hap' in s]
        files_hap_2 = [s for s in files if '_2.hap' in s]
        mal_parent_list = []
        mal_stem_list = []
        for f in files_hap_1:
            f_stem = ''.join(f.split('_')[:-1])
            matched = [s for s in files_hap_2 if f_stem + '_2.hap' in s]
            if len(matched) == 1:
                mal_parent_list.append('{}/{}'.format(self.mal_dir, f_stem))
                mal_stem_list.append(f_stem)

        throw = []
        #make seed of random numbers
        seed_txt = '{}__{}__{}'.format(self.fem_dir, self.mal_dir, self.seed)
        hash = hashlib.md5(seed_txt.encode("utf-8")).hexdigest()
        #6 digits are used for seed
        hash = int(hash[:6], 16)  
        random.seed(hash)
        for h in range(len(fem_parent_list)):
            p1_stem = fem_parent_list[h].split('/')[-1]
            for i in range(self.num):
                if self.num < 10:
                    num_str = '{:01}'.format(i+1)
                elif self.num < 100:
                    num_str = '{:02}'.format(i+1)
                else:
                    num_str = '{:03}'.format(i+1)
                stem = '{}x(rd)-{}'.format(p1_stem, num_str)
                rand = random.randrange(len(mal_parent_list))
                one = [self.map, fem_parent_list[h], mal_parent_list[rand],
                       stem, self.seed, self.out]
                throw.append(one)

        #Multi processing
        ##############################################
        with Pool(self.cpu) as p:
            p.map(recombination.recombination, throw)
        #############################################

def main():
    print(time_stamp(), 'Randomcross simulation started.', flush=True)

    prog = Randcross(args)
    prog.run()

    print(time_stamp(), 'Successfully finished.\n', flush=True)

if __name__ == '__main__':
    main()
