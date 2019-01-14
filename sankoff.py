# -*- coding: utf-8 -*-
"""
Created on Fri Dec 28 15:13:05 2018

@author: Sharon
"""

# A substitution cost matrix
#        A  C  G  T  Gap
cost = [[0, 4, 2, 4, 8], # A
        [4, 0, 4, 2, 8], # C
        [2, 4, 0, 4, 8], # G
        [4, 2, 4, 0, 8], # T
        [8, 8, 8, 8, 8]] #Gap

inf = float('Inf')  # Infinity
def join(node1, node2):
    """Join two nodes together using the global cost matrix."""
    this = [0, 0, 0, 0, 0]
    
    for i in range(0, 5):
        # Find min-cost change from node on the left
        min_left = inf
        for j in range(0, 5):
            this_cost = cost[i][j] + node1[j]
            min_left = min(min_left, this_cost)
        # Find min-cost change from node on the right
        min_right = inf
        for k in range(0, 5):
            this_cost = cost[i][k] + node2[k]
            min_right = min(min_right, this_cost)
        this[i] = min_left + min_right  
    return this
