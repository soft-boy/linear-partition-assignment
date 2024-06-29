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

phi = lambda x: x
plus = lambda x, y: x+y
mult = lambda x, y: x*y
defaultsemiring=(phi, plus, mult, 0, 1)

def linear_partition(rna, b=100, semiring=defaultsemiring, T=10000000):
  phi, plus, mult, zero, one = semiring

  n = len(rna)
  Q = [defaultdict(lambda: float(zero)) for _ in range(n+1)]
  for j in range(1, n+1): Q[j-1][j] = one

  for j in range(1, n+1):
    for i in Q[j-1]:
      Q[j][i] = plus(Q[j][i], mult(Q[j-1][i], phi(delta_eq)))
      
      pair = get_pair(rna, i-1, j)
      if pair in {"AU", "UA", "CG", "GC", "GU", "UG"}:
        for k in Q[i-2]:
          delta = phi(pow(e, -xi(rna, i-1, j)/(R*T)))
          delta = mult(delta, Q[j-1][i])
          delta = mult(delta, Q[i-2][k])
          Q[j][k] = plus(Q[j][k], delta)
    beam_prune(Q, j, b)

  return Q

def base_pairing_probs(rna, Q, T=10000000):
  n = len(rna)
  Q_hat = defaultdict(float)
  p = defaultdict(float)
  Q_hat[(1, n)] = 1

  for j in range(n, 0, -1):
    for i in Q[j-1]:
      Q_hat[(i, j-1)] += Q_hat[(i, j)]*delta_eq

      pair = get_pair(rna, i-1, j)
      if pair in {"AU", "UA", "CG", "GC", "GU", "UG"} and i-2 > -1:
        for k in Q[i-2]:
          energy = xi(rna, i-1, j)
          Q_hat[(k, i-2)] += Q_hat[(k, j)] * Q[j-1][i] * pow(e, -energy/(R*T))
          Q_hat[(i, j-1)] += Q_hat[(k, j)] * Q[i-2][k] * pow(e, -energy/(R*T))
          p[(i-1,j)] += (Q_hat[(k, j)] * Q[i-2][k] * math.exp(-energy/(R*T)) * Q[j-1][i])/Q[n][1]
          p[(j,i-1)] += (Q_hat[(k, j)] * Q[i-2][k] * math.exp(-energy/(R*T)) * Q[j-1][i])/Q[n][1]

  return p, Q_hat