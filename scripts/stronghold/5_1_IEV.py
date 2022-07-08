# 5_1_IEV
"""
See note
"""
import scripts.Rosalind as rs


### read data ###
with open("../../inputs/stronghold/rosalind_iev.txt", "r") as f:
    pops = list(map(int, f.readline().rstrip().split()))

### calculate dominant phenotype probability
# respective probability for dominant offspring
probs = [1, 1, 1, .75, .5, 0]
# number of offsprings per parent pair
offsprings = 2

res = sum([rs.binomial_sum(pr, offsprings, coef=True) * pop
           for pr, pop in zip(probs, pops)])

print(res)
