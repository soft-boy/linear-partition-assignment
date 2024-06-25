from collections import defaultdict
import math
from linear_partition import linear_partition

R = 8.31446
in_edges = lambda v: [] #TODO
w = lambda e: 0 #TODO
Z = lambda v: 0 #TODO

S = defaultdict(set)

def sum_edges(v, saving):
  s = 0
  for e in in_edges(v):
    _, subs = e
    s += math.exp((-w(e)/(R*T) * math.prod([Z(u) for u in subs])))
    if saving: S[v].add((e, s))
  
  return s

sample = lambda l: l[0] #TODO
f = lambda e: 0 #TODO

def sample_lazy(v, visited, Q):
  if v not in visited:
    sum_edges(v, True)
    visited.add(v)
  else:
    e = sample(S[v])
    _, subs = e
    #return f(e)(sample_lazy([(u, visited) for u in subs]))
  
def main_lazy(rna, k):
  n = len(rna)
  Q = linear_partition(rna)
  visited = set()
  for _ in range(k):
    sample_lazy((rna, 1, n), visited)