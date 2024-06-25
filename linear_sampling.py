from collections import defaultdict
import math
from linear_partition import linear_partition, xi

R = 8.31446
w = lambda e: 0 #TODO
Z = lambda v: 0 #TODO
pair_match = lambda x, i, j: xi(x, i, j) != 0

S = defaultdict(set)

def in_edges(v, Z):
  x, i, j = v
  if j == i-1:
    return ((x, i, i-1), [])
  if Z[j-1][i] != 0:
    return ((x, i, j), [(x, i, j-1)])
  
  for k_p in Z[j-1]:
    k = k_p-1
    if pair_match(x, k, j) and Z[k-1][i] != 0:
      return ((x, i, j), [(x, i, k-1), (x, k+1, j-1)])


def sum_edges(v, saving):
  s = 0
  for e in in_edges(v):
    _, subs = e
    s += math.exp((-w(e)/(R*T) * math.prod([Z(u) for u in subs])))
    if saving: S[v].add((e, s))
  
  return s

sample = lambda l: l[0] #TODO
f = lambda e: 0 #TODO

def sample_lazy(v, visited, Z):
  if v not in visited:
    sum_edges(v, True)
    visited.add(v)
  else:
    e = sample(S[v])
    _, subs = e
    #return f(e)(sample_lazy([(u, visited) for u in subs]))
  
def main_lazy(x, k):
  n = len(x)
  Z = linear_partition(x)
  visited = set()
  samples = []
  for _ in range(k):
    samples.append(sample_lazy((x, 1, n), visited, Z))

  return samples