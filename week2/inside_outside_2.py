from collections import defaultdict
from qselect import qselect
import math

e = 2.71828
R = 8.31446
delta_eq = 1

get_pair = lambda x, i, j: x[i-1]+x[j-1]

def xi(x, i, j):
  pair = get_pair(x, i, j)
  if pair in {"CG", "GC"}: return -3
  elif pair in {"AU", "UA"}: return -2
  elif pair in {"GU", "UG"}: return -1
  else: return 0

def top_b(candidates, b):
  if len(candidates) <= b: return candidates
  
  scores = [(score, i) for i, score in candidates.items()]
  b_score = qselect(b, scores)[0]

  return {k: v for k, v in candidates.items() if v >= b_score}

def beam_prune(Q, j, b):
  candidates = {}

  li = [i for i in Q[j]]
  for i in li:
    candidates[i] = Q[i-1][1]*Q[j][i]
  
  candidates = top_b(candidates, b)

  to_remove = [i for i in Q[j] if i not in candidates]
  for i in to_remove: del Q[j][i]

def inside_2(rna, b=100, T=10000000):
  n = len(rna)
  S_in = [defaultdict(lambda: 0) for _ in range(n+1)]
  P_in = [defaultdict(lambda: 0) for _ in range(n+1)]

  for j in range(1, n+1): S_in[j-1][j] = 1

  for j in range(1, n+1):
    for i in S_in[j-1]:
      S_in[j][i] += S_in[j-1][i] * delta_eq
      
      pair = get_pair(rna, i-1, j)
      if pair in {"AU", "UA", "CG", "GC", "GU", "UG"}:
        for k in S_in[i-2]:
          P_in[j][i-1] = S_in[j-1][i] * pow(e, -xi(rna, i-1, j)/(R*T))
          S_in[j][k] += P_in[j][i-1] * S_in[i-2][k]
    
    beam_prune(S_in, j, b)

  return S_in, P_in

def outside_2(rna, S_in, P_in, T=10000000):
  n = len(rna)
  S_out = defaultdict(float)
  P_out = defaultdict(float)
  p = defaultdict(float)
  S_out[(1, n)] = 1

  for j in range(n, 0, -1):
    for i in S_in[j-1]:
      S_out[(i, j-1)] += S_out[(i, j)] * delta_eq

      pair = get_pair(rna, i-1, j)
      if pair in {"AU", "UA", "CG", "GC", "GU", "UG"} and i-2 > -1:
        for k in S_in[i-2]:
          S_out[(k, i-2)] += S_out[(k, j)] * P_in[j][i-1]
          P_out[(i-1, j)] += S_out[(k, j)] * S_in[i-2][k]

          p[(i-1,j)] += P_out[(i-1, j)] * P_in[j][i-1] / S_in[n][1]
          p[(j,i-1)] = p[(i-1,j)]

  return p, S_out, P_out