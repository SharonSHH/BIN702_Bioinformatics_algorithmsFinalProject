# -*- coding: utf-8 -*-
"""
Created on Tue Dec 25 19:12:00 2018

@author: Sharon
"""
from Bio.Align.Applications import MuscleCommandline

in_file = "C:/BIN702/input.fasta"
out_file = "C:/BIN702/aligned.fasta"
muscle_exe = "muscle3.8.31_i86win32.exe"
try:
    cline = MuscleCommandline(muscle_exe, input=in_file, out=out_file)
    stdout, stderr = cline(in_file)
    from StringIO import StringIO
    from Bio import AlignIO
    align = AlignIO.read(StringIO(stdout), "fasta")
    AlignIO.close()
except ValueError:
    print "There is no valid number."
