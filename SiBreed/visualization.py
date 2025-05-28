import math
import sys

from matplotlib import patches, pyplot as plt
from SiBreed import rwhap
from SiBreed.utils import time_stamp

colors =['lightgray', 'black', 'red', 'blue', 'yellow', 
         'darkgreen', 'orange', 'pink', 'aqua', 'lightgreen',
         'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray',
         'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray',
         'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray', 'gray',
         'gray', 'gray']

def visualization(list):
    parent = list[0]
    char_1 = list[1]
    out = list[2]
    #Get length of chromosome
    chr = []
    chr_len = []
    hap = [rwhap.readhap('{}_1.hap'.format(parent)),
            rwhap.readhap('{}_2.hap'.format(parent))]
    format = hap[0][1]
    for i in range(len(hap[0]) - 2):
        chr.append(hap[0][i+2][0])
        seq = hap[0][i+2][1]
        chr_len.append(len(seq) * format)

    #number of digits in the length of the longest chromosome
    max_len = max(chr_len)
    digits = math.floor(math.log10(max_len))
    standard = 10**(digits)
    #if the longest chr length is 23098790, standard = 10000000
    if(max_len / standard < 2):
        standard = standard / 5
    elif(max_len / standard < 5):
        standard = int(standard / 2)
    #if the longest chr length is 23098790, standard = 5000000
    y_axis_at = range(0, standard*11, standard)
    y_axis_lab = []
    if(standard >= 100000):
        st_lab = standard/1000000
        sign = 'M'
    elif(standard >= 100):
        st_lab = standard/1000
        sign = 'K'
    else:
        st_lab = standard
        sign = 'bp'
    for i in range(11):
        y_axis_lab.append('{}{}'.format(round(st_lab * i, 1), sign))

    #Make figure
    hap = [rwhap.readhap('{}_1.hap'.format(parent)),
            rwhap.readhap('{}_2.hap'.format(parent))]
    chr_num = len(hap[0]) - 2
    # Create a figure
    fig = plt.figure(figsize=(5, 5), dpi=144)
    ax = fig.add_subplot(111,
                        xlim=[-1, len(chr)],
                        xticks=range(len(chr)),
                        xticklabels=chr,
                        xlabel="Chromosome",
                        ylim=[max_len*1.05, -max_len*0.05],
                        yticks=y_axis_at,
                        yticklabels=y_axis_lab,
                        ylabel="Position")
    plt.subplots_adjust(left=0.15, right=0.95, bottom=0.15, top=0.95)
    plt.xticks(rotation=45)
    plt.xlim(-1, len(chr))
    plt.ylim(max_len*1.05, -max_len*0.05)

    legends = []
    for i in range(len(char_1)):
        legends.append(patches.Patch(color=colors[i], 
                                        label=char_1[i]))
    ax.legend(handles=legends, loc="lower right")

    #each chromosome
    for i in range(chr_num):
        adj = [-0.4, 0]
        for j in range(2):
            str = hap[j][i + 2][1]  #sequence
            current_char_index = 0
            top = 0
            for k in range(len(str)):
                if str[k] != char_1[current_char_index]:
                    r = patches.Rectangle(xy=(i + adj[j], top), width=0.4, 
                        height=k * format - top, ec=None, 
                        fc=colors[current_char_index], fill=True)
                    ax.add_patch(r)
                    flag = True
                    for l in range(len(char_1)):
                        if str[k] == char_1[l]:
                            current_char_index = l
                            top = k * format
                            flag = False
                            break
                    if flag:
                        print(time_stamp(),
                              'Haploid files have unrecognized characterr.',
                              flush=True)
                        sys.exit(1) 
            #To the terminal
            r = patches.Rectangle(xy=(i + adj[j], top), width=0.4, 
                        height=chr_len[i] - top, ec=None, 
                        fc=colors[current_char_index], fill=True)
            ax.add_patch(r)

            #Draw rectangle of chromosome
            r = patches.Rectangle(xy=(i + adj[j], 0), width=0.4,
                height=chr_len[i], ec='black', fill=False)
            ax.add_patch(r)

    # Save figure
    file = '{}/{}.png'.format(out, parent.split('/')[-1])
    fig.savefig(file, dpi=144)

    # Release memory
    plt.clf()
    plt.close()
