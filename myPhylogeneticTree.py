# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 19:12:00 2018

@author: Sharon
"""

import sankoff
from Bio import SeqIO
import os
import copy
import re


def getMultiAlignment():
    pathToFile = os.path.join("aligned.fasta")
    names = []
    dict = {}
    for seq_record in SeqIO.parse(pathToFile, """fasta"""):
        names.append(seq_record.name)
        dict[seq_record.name] = seq_record.seq
    return dict, names


def getSankoffValue(val):
        inf = float('Inf') 
        if val.upper() == 'A':
            return [0, inf, inf, inf, inf]
        elif val.upper() == 'C':
            return [inf, 0, inf, inf, inf]
        elif val.upper() == 'G':
            return [inf, inf, 0, inf, inf]
        elif val.upper() == 'T':
            return [inf, inf, inf, 0, inf]
        elif val.upper() == '-':
            return [inf, inf, inf, inf, 0]
        
        
def selectMin(fiveElements):
        min_value = min(fiveElements)
        element = fiveElements.index(min_value)
        if element == 0:
            return 'A'
        elif element == 1:
            return 'C'
        elif element == 2:
            return 'G'
        elif element == 3:
            return 'T'
        else:
            return '-'
        

#Build the trees with only two nodes, left node and right node            
def createTree(list, dict, index):
        sub_dict = {}
        nodeList = [] #save the root of trees
        stack = []
        count = 0
        
        while True:
            #There is no data in the list, it means trees have already built
            if len(list)==0:
                return nodeList, sub_dict, count
            # input to stack when '(' and strings
            if list[0] == '(':
                stack.append(list[0])
                list.pop(0)
            elif list[0] != ')': # String
                if stack[-1] == '(': #The first element is (, input list[0]
                    stack.append(list[0])
                    list.pop(0)
                else:
                    lvalue = getSankoffValue(dict[stack[-1]][index])
                    rvalue = getSankoffValue(dict[list[0]][index])
                    rootValue = sankoff.join(lvalue, rvalue)
                    rootName = 'root' + str(count) 
                    sub_dict[rootName] = rootValue
                    count += 1
                    nodeList.append(rootName)
                    preRoot = rootName
                    list.pop(0)
                    stack.pop(-1) #pop the last element
            else: # input is ')'
                if stack[-1] == '(': #The newest element is (, popup one element from list and stack 
                    stack.pop(-1)
                    list.pop(0)
                else: #The newest element of stack is string, build new node and popup two elements of stack
                    lvalue = sub_dict[preRoot]
                    nodeList.pop(-1)
                    rvalue = getSankoffValue(dict[stack[-1]][index])
                    rootValue = sankoff.join(lvalue, rvalue)
                    rootName = 'root' + str(count)
                    sub_dict[rootName] = rootValue
                    nodeList.append(rootName)
                    count += 1
                    preRoot = rootName
                    list.pop(0)
                    stack.pop(-1)
                    stack.pop(-1)
        

#Scan the roots in nodeList and combine them into one binary tree    
def binaryTree(nodeList, sub_dict, count):
        newPreRoot = ''
        if len(nodeList) == 1:
            fiveElements = sub_dict[nodeList[0]]
            ancestralSeq = selectMin(fiveElements)
            nodeList.pop(0)
            return ancestralSeq
        for i in range(len(nodeList) - 1):
            if i == 0:
                lvalue = sub_dict[nodeList[0]]
                rvalue = sub_dict[nodeList[1]] 
            else:
                lvalue = sub_dict[newPreRoot]
                rvalue = sub_dict(nodeList[i+1])
        rootValue = sankoff.join(lvalue, rvalue)
        rootName = 'root' + str(count)
        sub_dict[rootName] = rootValue
        newPreRoot = rootName
        count += 1
        fiveElements = sub_dict[newPreRoot]
        ancestralSeq = selectMin(fiveElements)
        return ancestralSeq


if __name__ == '__main__':
    f = open("C:/BIN702/phylogeny.txt_phyml_tree.txt" , "r")
    treedata = f.readlines()
    pattern = re.compile(r'\(|\)|[a-zA-Z]+')
    Original = pattern.findall(treedata[0])
    print(Original)

    dict, names = getMultiAlignment()
    ancestralSeq = ''
     
    for i in range(len(dict[names[0]])):
        list = copy.copy(Original)
        nodeList, sub_dict, count = createTree(list, dict, i)
        ancestralSeq += binaryTree(nodeList, sub_dict, count) 

    if (len(ancestralSeq.replace("-", "")) != len(ancestralSeq)):
        ancestralSeq += '\nAncestral Sequence without gap:\n' + ancestralSeq.replace("-", "")

    with open("C:/BIN702/Ancestral Sequence.txt", "w") as text_file:
        text_file.write("Ancestral Sequence:\n{0}".format(ancestralSeq.strip()))
    
    print('\n')
    print(ancestralSeq)