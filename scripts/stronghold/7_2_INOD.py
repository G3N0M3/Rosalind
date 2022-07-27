# 7_2_INOD

### read data ###
with open("../../inputs/stronghold/rosalind_inod.txt", "r") as f:
    n = int(f.readline().rstrip())

### calculate ###
"""
When n = 3, internal nodes that satisfy conditions are 1
When n increase by 1, number of nodes also increase by one.
Thus, when number of nodes that satisfy conditions is "y"
y = n - 2 
"""

print(n - 2)
