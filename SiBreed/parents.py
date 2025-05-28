#!/usr/bin/env python3
import math
import os
from SiBreed import rwhap
from SiBreed.params import Params
from SiBreed.utils import command, time_stamp

pm = Params('parents')
args = pm.set_options()

symbols = ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
           'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',
           'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',
           'u', 'v', 'w', 'y', 'z']

class Parents(object):
    def __init__(self, args):
        self.args = args
        self.out = args.out
        self.num = args.num
        self.bp = args.base_per_character

        #Check args
        pm.parents_check_args(args)

        #Read map file
        self.map = rwhap.readmap(args.map)

        #Make output directory
        os.mkdir('{}'.format(self.out))

        #Write command infomation
        command(self.out)

    def run(self):
        for h in range(1,3):
            for i in range(self.num):
                data = []
                data.append('P{}_{}'.format(symbols[i], h))
                data.append(self.bp)
                chrmo = []
                current_chr = ''
                for j in range(1, len(self.map)): #map[0] is header
                    if current_chr == '':
                        current_chr = self.map[j][0]
                    elif current_chr != self.map[j][0]:
                        chrmo.append(self.map[j-1][0])
                        char_len = math.ceil(int(self.map[j-1][1]) / self.bp)
                        chrmo.append(symbols[i] * char_len)
                        data.append(chrmo)

                        chrmo = []
                        current_chr = self.map[j][0]
                        
                chrmo.append(self.map[j][0])
                char_len = math.ceil(int(self.map[j][1]) / self.bp)
                chrmo.append(symbols[i] * char_len)
                data.append(chrmo)

                rwhap.writehap(data, '{}/P{}_{}.hap'.format(self.out, symbols[i], h))

def main():
    print(time_stamp(), 'Generating haploid files of parents started.', flush=True)

    prog = Parents(args)
    prog.run()

    print(time_stamp(), 'Successfully finished.\n', flush=True)

if __name__ == '__main__':
    main()
