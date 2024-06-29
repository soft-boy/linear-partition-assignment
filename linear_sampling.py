from collections import defaultdict
import math
import random
from linear_partition import linear_partition, base_pairing_probs, xi, get_pair

R = 8.31446
T = 1
S = defaultdict(list)

def in_edges(v, Z):
  x, i, j = v
  edges= []

  if j == i-1:
    return [((x, i, i-1), [])]
  if Z[j-1][i] != 0:
    edges.append(((x, i, j), [(x, i, j-1)]))
  
  for k_p in Z[j-1]:
    k = k_p-1
    pair = get_pair(x, k, j)
    if pair in {"AU", "UA", "CG", "GC", "GU", "UG"} and Z[k-1][i] != 0:
      edges.append(((x, i, j), [(x, i, k-1), (x, k+1, j-1)]))
  
  return edges

def sum_edges(v, saving, Z):
  s = 0
  for e in in_edges(v, Z):
    v, subs = e
    x, i, j = v
    s += math.exp((-xi(x, i, j)/(R*T) * math.prod([Z[j][i] for (_, i, j) in subs])))
    if saving: S[v].append((e, s))
  
  return s

def binary_search(lst, p):
  if not lst:
    return None

  low = 0
  high = len(lst) - 1
  result = None

  while low <= high:
    mid = (low + high) // 2
    result = lst[mid]  
    if lst[mid][1] > p:
      high = mid - 1 
      if high >= 0 and lst[high][1] < p:
        return result[0]  
    else:
      low = mid + 1 

  return result[0]  

def sample_lazy(v, visited, Z):
  _, i, j = v

  if v not in visited:
    sum_edges(v, True, Z)
    visited.add(v)

  h_edges= S[v]
    
  if h_edges:
    p = random.uniform(0, Z[j][i])
    e = binary_search(h_edges, p)
    _, subs = e
      
    if len(subs)==1:
      return sample_lazy(subs[0], visited, Z)+'.'
    elif len(subs)==2:
      return sample_lazy(subs[0], visited, Z)+'('+sample_lazy(subs[1], visited, Z)+ ')'
    else:
      return ''
  else:
    return ''
  
def parse_structure(structure):
  stack = []
  pairs = {}
  for index, char in enumerate(structure):
    if char == '(':
      stack.append(index)
    elif char == ')':
      if stack:
        opening_index = stack.pop()
        pairs[opening_index] = index
  return pairs
  
def aggregate_pairs(structures):
  pair_counts = {}
  for structure in structures:
    pairs = parse_structure(structure)
    for i, j in pairs.items():
      if (i, j) not in pair_counts:
        pair_counts[(i, j)] = 0
      pair_counts[(i, j)] += 1
  return pair_counts

def calculate_probabilities(pair_counts, total_structures):
  probabilities = defaultdict(float)
  for pair, count in pair_counts.items():
    probabilities[pair[0]+1, pair[1]+1] = float(count) / total_structures
    probabilities[pair[1]+1, pair[0]+1] = float(count) / total_structures
  return probabilities
  
def test(x, k):
  n = len(x)
  Z = linear_partition(x, T=T)
  visited = set()
  samples = []
  for _ in range(k):
    samples.append(sample_lazy((x, 1, n), visited, Z))

  
  pair_count= aggregate_pairs(samples)
  p, _ = base_pairing_probs(x, Z, T=T)
  sample_p= calculate_probabilities(pair_count, len(samples))

  print("The probability from the parititon function:")
  print()
  print(dict(p))

  print()
  print("The probability from sampling:")
  print()
  print(sample_p)

# test('AGACU', 100000)