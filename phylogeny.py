# -*- coding: utf-8 -*-
"""
Created on Fri Dec 28 12:32:33 2018

@author: Sharon
"""

from Bio.Phylo.Applications import PhymlCommandline
from Bio import Phylo
import os


def insert(original, new, pos):
    return original[:pos] + new + original[pos:]


def fastaToPhylogeny():  
    seq={}
    count = 0
    num = 0
    src = 'C:/BIN702/input.fasta'
    dst = 'C:/BIN702/phylogeny.txt' 
    f = open(src, 'r')
    content = ''
    for line in f:
        if line.startswith('>'):
            name=line.replace('>','').split()[0]
            seq[name]=''
            content += name + '\n'
        else:
            seq[name]+=line.replace('\n','').strip()
            if count == 1:
                num = len(seq[name])
            content += seq[name] + '\n' 
        count = len(seq)
    f.close()
    newstr = str(count) + ' ' + str(num) + '\n'
    content = insert(content, newstr, 0)
    fwrite = open(dst, 'w')
    fwrite.write(content)
    fwrite.close()
    
    
if __name__ == '__main__':
    fastaToPhylogeny()
    EX_PHYLIP = "C:/BIN702/phylogeny.txt"
    out_file = "C:/BIN702/phyml_tree.fasta"
    phyml_exe = "PhyML-3.1_win32.exe"
    os.environ['LANG'] = 'C'
    cmd = PhymlCommandline(phyml_exe, input=EX_PHYLIP, datatype='aa')
    
    try:
        out, err = cmd()
        # Check the output tree
        outfname = EX_PHYLIP + '_phyml_tree.txt'
        if not os.path.isfile(outfname):
            # NB: Briefly, PhyML dropped the .txt suffix (#919)
            outfname = outfname[:-4]
        tree = Phylo.read(outfname, 'newick')
        
    except Exception as exc:
        print("PhyML wrapper error: %s" % exc)
