#!/usr/bin/env python3
import os
import sys
from SiBreed import rwhap
from SiBreed.params import Params
from SiBreed.utils import time_stamp

pm = Params('genostat')
args = pm.set_options()

class GenoStat(object):
    def __init__(self, args):
        self.args = args
        self.dir = args.directory

        #Check args
        pm.genostat_check_args(args)

        self.char_1 = args.character.split(',')
        self.char_2 = []
        for i in range(len(self.char_1)):
            for j in range(len(self.char_1)):
                self.char_2.append('{}/{}'.format(self.char_1[i],
                                                  self.char_1[j]))

        #Make output directory
        self.out = '{}_stat'.format(self.dir)
        os.mkdir(self.out)

    def run(self):
        files = os.listdir(self.dir)
        files_hap_1 = [s for s in files if '_1.hap' in s]
        files_hap_2 = [s for s in files if '_2.hap' in s]
        parent_list = []
        stem_list = []
        sample_num = 0
        for f in files_hap_1:
            f_stem = '_'.join(f.split('_')[:-1])
            matched = [s for s in files_hap_2 if f_stem + '_2.hap' in s]
            if len(matched) == 1:
                parent_list.append('{}/{}'.format(self.dir, f_stem))
                stem_list.append(f_stem)
                sample_num = sample_num + 1
        if sample_num < 1:
            print(time_stamp(),
                  'Disgnated directory does not cantain valid haploid files.',
                  flush=True)
            sys.exit(1) 

        n = len(self.char_1)
        sum_1 = n * [0]
        sum_2 = (n * n) * [0]

        for h in range(len(parent_list)):
            hap = [rwhap.readhap('{}_1.hap'.format(parent_list[h])),
                   rwhap.readhap('{}_2.hap'.format(parent_list[h]))]
            chr_num = len(hap[0]) - 2
            list_a = [hap[0][0]]
            cnt_all_1 = n * [0]
            cnt_all_2 = (n * n) * [0]

            #stats of each chromosome
            for i in range(chr_num):
                list_b = [hap[0][i + 2][0]] #[0]:chromosome name
                str_1 = hap[0][i + 2][1]
                str_2 = hap[1][i + 2][1]
                cnt_1 = n * [0]
                cnt_2 = (n * n) * [0]
                for j in range(len(str_1)):
                    x = -1
                    y = -1
                    for k in range(n):
                        if str_1[j] == self.char_1[k]:
                            x = k
                            break
                    for k in range(n):
                        if str_2[j] == self.char_1[k]:
                            y = k
                            break
                    if x == -1 or y == -1:
                        print(time_stamp(),
                              'Haploid files have unrecognized characterr.',
                              flush=True)
                        sys.exit(1) 
                    cnt_1[x] = cnt_1[x] + 1
                    cnt_1[y] = cnt_1[y] + 1
                    cnt_2[x * n + y] = cnt_2[x * n + y] + 1
                list_b.append(self.char_1 + self.char_2) #[1]:genotype name
                list_b.append(cnt_1 + cnt_2) #[2]:count of genotype
                par_1 = [c / sum(cnt_1) for c in cnt_1]
                par_2 = [c / sum(cnt_2) for c in cnt_2]
                list_b.append(par_1 + par_2) #[3]:percentage of genotype
                list_a.append(list_b)
                cnt_all_1 = [ca + c for ca, c in zip(cnt_all_1, cnt_1)]
                cnt_all_2 = [ca + c for ca, c in zip(cnt_all_2, cnt_2)]

            #stats of whole genome
            list_b_all = ['All'] #[0]:chromosome name
            list_b_all.append(self.char_1 + self.char_2) #[1]:genotype name
            list_b_all.append(cnt_all_1 + cnt_all_2) #[2]:count of genotype
            par_all_1 = [ca / sum(cnt_all_1) for ca in cnt_all_1]
            par_all_2 = [ca / sum(cnt_all_2) for ca in cnt_all_2]
            list_b_all.append(par_all_1 + par_all_2) #[3]:percentage of genotype
            list_a.append(list_b_all)

            rwhap.writestat(list_a, '{}/{}.stat'.format(self.out, stem_list[h]))

            #For summary of directory
            sum_1 = [cs + ca for cs, ca in zip(sum_1, cnt_all_1)]
            sum_2 = [cs + ca for cs, ca in zip(sum_2, cnt_all_2)]

        #summary of directory
        list_s = [self.dir.split('/')[-1]]
        list_s_sub = ['Summary'] #[0]:chromosome name
        list_s_sub.append(self.char_1 + self.char_2) #[1]:genotype name
        list_s_sub.append(sum_1 + sum_2) #[2]:count of genotype
        par_sum_1 = [cs / sum(sum_1) for cs in sum_1]
        par_sum_2 = [cs / sum(sum_2) for cs in sum_2]
        list_s_sub.append(par_sum_1 + par_sum_2) #[3]:percentage of genotype
        list_s.append(list_s_sub)

        rwhap.writestat(list_s, 
                        '{}/{}_summary.stat'.format(self.out, 
                                                    self.dir.split('/')[-1]))


def main():
    print(time_stamp(), 'Genotype stat started.', flush=True)

    prog = GenoStat(args)
    prog.run()

    print(time_stamp(), 'Successfully finished.\n', flush=True)

if __name__ == '__main__':
    main()
