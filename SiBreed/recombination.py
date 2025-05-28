import hashlib
import math
import random
import sys
from SiBreed import rwhap
from SiBreed.utils import time_stamp

def recombination(list):
    #Map file
    map = list[0]
    #Read haploid files
    hap = [[rwhap.readhap('{}_1.hap'.format(list[1])),
            rwhap.readhap('{}_2.hap'.format(list[1]))],
           [rwhap.readhap('{}_1.hap'.format(list[2])),
            rwhap.readhap('{}_2.hap'.format(list[2]))]]
    #check format
    if not(hap[0][0][1] == hap[0][1][1] == 
           hap[1][0][1] == hap[1][1][1]):
        print(time_stamp(),
            'Format of haploid files are not same among parents.',
            flush=True)
        sys.exit(1) 
    else:
        format = int(hap[0][0][1]) 
    
    #Output file name
    stem = list[3]
    #Seed
    seed = list[4]
    #Output directory
    outdir = list[5]
            
    for h in range(1,3):
        #make seed of random numbers
        seed_txt = '{}__{}__{}__{}'.format(outdir, h, stem, seed)
        hash = hashlib.md5(seed_txt.encode("utf-8")).hexdigest()
        #6 digits are used for seed
        hash = int(hash[:6], 16)  
        random.seed(hash)

        #Simulate recombinanation events
        events = []
        current_chr = ''
        map_sub = []
        for j in range(1, len(map)): #map[0] is header
            if current_chr != map[j][0]:
                current_chr = map[j][0]
                if len(map_sub) != 0:
                    events = events + simrecombi(map_sub)
                    map_sub = []
            map_sub.append(map[j])
        events = events + simrecombi(map_sub)

        #Apply recombination events
        new_gamete = []
        new_gamete.append('{}_{}'.format(stem, h)) #This is the name
        new_gamete.append(format)
        for i in range(2, len(hap[h-1][0])):
            current_chr = hap[h-1][0][i][0]
            current_seq = [hap[h-1][0][i][1], hap[h-1][1][i][1]]
            events_sub = []
            new_seq = ''
            for j in range(len(events)):
                if events[j][0] == current_chr:
                    events_sub.append(events[j])
            for j in range(len(events_sub)):
                if j == 0:
                    new_seq = current_seq[events_sub[j][2]]
                else:
                    char_pos = math.floor(events_sub[j][1] / format) 
                    new_seq = new_seq[:char_pos] + current_seq[events_sub[j][2]][char_pos:]

            new_gamete.append([current_chr, new_seq])

        rwhap.writehap(new_gamete, '{}/{}_{}.hap'.format(outdir, stem, h))

def simrecombi(map):
    cM_start = float(map[0][2])
    cM_end = float(map[len(map) - 1][2])
    cM_dif = cM_end - cM_start

    #Recombination events occurs based on Poisson distribution
    rand = random.random()
    c = 0
    prob = math.e ** (-2 * cM_dif / 100)
    while(rand > prob and c < 20): #limit of recombination: 20
        c = c + 1
        prob = prob + ((2 * cM_dif / 100) ** c) * (math.e ** (-2 * cM_dif / 100)) / math.factorial(c)
   
    current_hap = int(round(random.random()))
    events = [[map[0][0], 1, current_hap]]
    if c > 0:
        cross_pos = []
        for i in range(c):
            cross_pos.append(random.uniform(cM_start, cM_end))
        cross_pos.sort()
        for pos in cross_pos:
            for j in range(len(map) - 1):
                if float(map[j][2]) <= pos and float(map[j+1][2]) >= pos:
                    slope = (int(map[j+1][1]) - int(map[j][1])) / (float(map[j+1][2]) - float(map[j][2]))
                    cross_pysical = int(round(slope * (float(map[j+1][2]) - pos))) + int(map[j][1])
                    #if current_hap == 0:
                    #    current_hap = 1
                    #else:
                    #    current_hap = 0
                    current_hap = int(round(random.random()))
                    events.append([map[j][0], cross_pysical, current_hap])
                    break
    return events