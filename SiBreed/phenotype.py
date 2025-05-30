#!/usr/bin/env python3
import csv
import hashlib
from multiprocessing import Pool
import os
import random
import sys
from SiBreed import rwhap
from SiBreed.params import Params
from SiBreed.utils import command, time_stamp

pm = Params('phenotype')
args = pm.set_options()

class Phenotype(object):
    def __init__(self, args):
        self.args = args
        self.qtl = args.qtl
        self.dir = args.directory
        self.sub = args.sub
        self.center = args.center
        self.err = args.error
        self.seed = args.seed
        self.cpu = args.cpu

        #Check args
        pm.phenotype_check_args(args)

        #Make output directory
        self.out = '{}_{}'.format(self.dir, self.sub)
        os.mkdir(self.out)

        #Write command infomation
        command(self.out)

    def run(self):
        throw_list = [] #for paralel processing

        #read qtl file
        qtl_table = rwhap.readqtl(self.qtl)

        #search valid haploid sets
        files = os.listdir(self.dir)
        files_hap_1 = [s for s in files if '_1.hap' in s]
        files_hap_2 = [s for s in files if '_2.hap' in s]  
        sample_num = 0
        for f in files_hap_1:
            f_stem = '_'.join(f.split('_')[:-1])
            matched = [s for s in files_hap_2 if f_stem + '_2.hap' in s]
            if len(matched) == 1:
                throw_list.append(['{}/{}'.format(self.dir, f_stem), qtl_table])
                sample_num = sample_num + 1
        if sample_num < 1:
            print(time_stamp(),
                  'Disgnated directory does not cantain valid haploid files.',
                  flush=True)
            sys.exit(1) 

        #Multi processing
        ##############################################
        with Pool(self.cpu) as p:
            result = p.map(calcpheno, throw_list)
        #############################################

        result = sorted(result, reverse=False, key=lambda x: x[0])

        #make seed of random numbers
        seed_txt = '{}__{}__{}'.format(self.dir, self.sub, self.seed)
        hash = hashlib.md5(seed_txt.encode("utf-8")).hexdigest()
        #6 digits are used for seed
        hash = int(hash[:6], 16)  
        random.seed(hash)

        for i in range(len(result)):
            result[i][1] = result[i][1] + self.center
            if self.err > 0:
                rand = random.gauss(mu=0.0, sigma=self.err)
                result[i][1] = result[i][1] + rand

        with open('{}/result_{}.tsv'.format(self.out, self.sub), 'w') as f:
            writer = csv.writer(f, delimiter='\t')
            writer.writerows(result)



def calcpheno(list):
    source = list[0]
    qtl = list[1]

    hap = [rwhap.readhap('{}_1.hap'.format(source)),
            rwhap.readhap('{}_2.hap'.format(source))]
    name = '_'.join(hap[0][0].split(sep='_')[:-1])
    format = int(hap[0][1])

    #This is base of phenotypic value
    value = 0.0

    current_chr = ''
    current_index = 0
    for i in range(1, len(qtl)): #first row is header
        if current_chr != qtl[i][0]:
            for j in range(2, len(hap[0])):
                if qtl[i][0] == hap[0][j][0]:
                    current_chr = qtl[i][0]
                    current_index = j
        alt = qtl[i][4].split(sep=',')
        pos = int(qtl[i][1]) // format
        char_1 = hap[0][current_index][1][pos]
        char_2 = hap[1][current_index][1][pos]
        flag_1 = False
        flag_2 = False
        for a in alt:
            if a == char_1:
                flag_1 = True
            if a == char_2:
                flag_2 = True
        if flag_1 and flag_2:
            value = value + float(qtl[i][2]) * 2
        elif flag_1 or flag_2:
            value = value + float(qtl[i][2]) + float(qtl[i][3])
        else:
            pass

    return [name, value]

def main():
    print(time_stamp(), 'Phenotyping started.', flush=True)

    prog = Phenotype(args)
    prog.run()

    print(time_stamp(), 'Successfully finished.\n', flush=True)

if __name__ == '__main__':
    main()
