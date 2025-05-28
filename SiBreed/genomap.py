#!/usr/bin/env python3
import csv
import os
import sys
from SiBreed import rwhap
from SiBreed.params import Params
from SiBreed.utils import time_stamp

pm = Params('genomap')
args = pm.set_options()

class GenoMap(object):
    def __init__(self, args):
        self.args = args
        self.dir = args.directory
        self.num = args.num

        #Check args
        pm.genomap_check_args(args)

        self.char_1 = args.character.split(',')
        if len(self.char_1) != 2:
            print(time_stamp(),
                  'It must be bi-parent population.',
                  flush=True)
            sys.exit(1) 

        #Make output directory
        self.out = '{}_map'.format(self.dir)
        os.mkdir(self.out)

    def run(self):
        files = os.listdir(self.dir)
        files_hap_1 = [s for s in files if '_1.hap' in s]
        files_hap_2 = [s for s in files if '_2.hap' in s]
        parent_list = []
        sample_num = 0
        for f in files_hap_1:
            f_stem = ''.join(f.split('_')[:-1])
            matched = [s for s in files_hap_2 if f_stem + '_2.hap' in s]
            if len(matched) == 1:
                parent_list.append('{}/{}'.format(self.dir, f_stem))
                sample_num = sample_num + 1
        if sample_num < 1:
            print(time_stamp(),
                  'Disgnated directory does not cantain valid haploid files.',
                  flush=True)
            sys.exit(1) 
        parent_list.sort()

        #row 1 and row 2
        header_1 = ['id']
        header_2 = ['']
        hap = [rwhap.readhap('{}_1.hap'.format(parent_list[0])),
               rwhap.readhap('{}_2.hap'.format(parent_list[0]))]
        for i in range(len(hap[0]) - 2):
            chr = hap[0][i+2][0]
            seq = hap[0][i+2][1]
            marker_num = (len(seq) // self.num) + 1
            for j in range(marker_num):
                header_1.append('{}_{:04}'.format(chr, j + 1))
                header_2.append(i + 1)

        rows = [header_1, header_2]
        for h in range(len(parent_list)):
            hap = [rwhap.readhap('{}_1.hap'.format(parent_list[h])),
                   rwhap.readhap('{}_2.hap'.format(parent_list[h]))]
            chr_num = len(hap[0]) - 2
            row = [hap[0][0]]

            #each chromosome
            for i in range(chr_num):
                str_1 = hap[0][i + 2][1]
                str_2 = hap[1][i + 2][1]

                for j in range(len(str_1)):
                    if j % self.num == 0:
                        if (str_1[j] == self.char_1[0] and 
                            str_2[j] == self.char_1[0]):
                            row.append('A')
                        elif (str_1[j] == self.char_1[0] and 
                              str_2[j] == self.char_1[1]):
                            row.append('H')
                        elif (str_1[j] == self.char_1[1] and 
                              str_2[j] == self.char_1[0]):
                            row.append('H')
                        elif (str_1[j] == self.char_1[1] and 
                              str_2[j] == self.char_1[1]):
                            row.append('B')
                        else:
                            print(time_stamp(),
                              'Haploid files have unrecognized characterr.',
                              flush=True)
                            sys.exit(1) 
            rows.append(row)

        with open('{}/{}_map.csv'.format(self.out, self.dir.split('/')[-1]),
                  'w') as f:
            writer = csv.writer(f)
            writer.writerows(rows)


def main():
    print(time_stamp(), 'Making genotype map started.', flush=True)

    prog = GenoMap(args)
    prog.run()

    print(time_stamp(), 'Successfully finished.\n', flush=True)

if __name__ == '__main__':
    main()
