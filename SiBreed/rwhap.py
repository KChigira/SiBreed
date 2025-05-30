import csv
import sys
from SiBreed.utils import time_stamp

def readmap(filename):
    #Read map file
    with open(filename, 'r') as f:
        map = list(csv.reader(f, delimiter='\t'))

    #Check format
    each_len = [len(v) for v in map]
    if max(each_len) != 3 or min(each_len) != 3:
        print(time_stamp(),
            'Input map file is malformatted.',
            flush=True)
        sys.exit(1) 
    return map

def readqtl(filename):
    #Read map file
    with open(filename, 'r') as f:
        map = list(csv.reader(f, delimiter='\t'))

    #Check format
    each_len = [len(v) for v in map]
    if max(each_len) != 5 or min(each_len) != 5:
        print(time_stamp(),
            'Input map file is malformatted.',
            flush=True)
        sys.exit(1) 
    return map

def readhap(filename):
    output = []
    chrom = []
    seq = ''
    with open(filename, 'r') as f:
        #Line 1
        l = f.readline()
        l = l.rstrip()
        if len(l) <= 6: error()
        if l[:6] == '#name=':
            output.append(l[6:])
        else:
            error()
        
        #Line 2
        l = f.readline()
        l = l.rstrip()
        if len(l) <= 8: error()
        if l[:8] == '#format=':
            try:
                output.append(int(l[8:]))
            except ValueError:
                error()
        else:
            error()

        #Line 3~
        lines = f.readlines()
        for l in lines:
            l = l.rstrip()
            if len(l) < 1: error()
            if l[:1] == '>':
                if seq != '':
                    chrom.append(seq)
                    if len(chrom) == 2:
                        output.append(chrom)
                        chrom = []
                    else:
                        error()
                chrom.append(l[1:])
                seq = ''
            else:
                seq = ''.join([seq, l])
        chrom.append(seq)
        if len(chrom) == 2:
            output.append(chrom)
        else:
            error()

    return output

def writehap(object, filename):
    with open(filename, 'w') as f:
        #Line 1
        if len(object) < 3: error()
        name = object[0]
        if isinstance(name, str):
            f.write('#name={}\n'.format(name)) 
        else:
            error()

        #Line 2
        format = object[1]
        if isinstance(format, int):
            f.write('#format={}\n'.format(format))
        else:
            error()

        #Line 3~
        for i in range(2, len(object)):
            if not(isinstance(object[i], list) and len(object[i]) == 2):
                error()
            if isinstance(object[i][0], str):
                f.write('>{}\n'.format(object[i][0])) 
            else:
                error()
            if isinstance(object[i][1], str):
                for j in range(len(object[i][1]) // 60):
                    f.write('{}\n'.format(object[i][1][(60 * j):(60 * (j+1))]))
                if len(object[i][1]) % 60 != 0:
                    f.write('{}\n'.format(object[i][1][(60 * (j+1)):]))
            else:
                error()

def writestat(object, filename):
    with open(filename, 'w') as f:
        #Line 1
        if len(object) < 2: error()
        name = object[0]
        if isinstance(name, str):
            f.write('#name={}\n'.format(name)) 
        else:
            error()

        #Line 2~
        for i in range(1, len(object)):                
            if isinstance(object[i][0], str):
                f.write('>{}\n'.format(object[i][0])) 
                f.write('genotype\tcount\tpercentage\n') 
            else:
                error()

            if len(object[i]) == 4:
                for j in range(len(object[i][1])):
                    f.write('{}\t{}\t{}\n'.format(object[i][1][j], 
                                                  object[i][2][j], 
                                                  object[i][3][j]))
            else:
                error()

def error():
    print(time_stamp(), 'Haploid file or object is malformatted.', flush=True)
    sys.exit(1) 
