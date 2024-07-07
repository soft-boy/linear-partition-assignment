from linear_partition import *
from linear_sampling import *
import numpy as np
np.set_printoptions(precision=2)

T = 1
k = 100000

def test(rna):
  print('RNA:', rna)
  print()

  n = len(rna)
  Q = linear_partition(rna, T=T)
  p, _ = base_pairing_probs(rna, Q, T=T)
  
  visited = set()
  samples = []
  for _ in range(k):
    samples.append(sample_lazy((rna, 1, n), visited, Q))

  pair_count= aggregate_pairs(samples)
  sample_p, _ = base_pairing_probs(rna, Q, T=T)
  sample_p = calculate_probabilities(pair_count, len(samples))

  p_arr = np.zeros((n, n))
  sample_p_arr = np.zeros((n, n))

  for i in range(1, n+1):
    for j in range(1, n+1):
      prob = p[(i, j)]
      sample_prob = sample_p[(i, j)]

      if prob < 0.001:
        prob = 0.0
      if sample_prob < 0.001:
        sample_prob = 0.0

      p_arr[i-1][j-1] = prob
      p_arr[j-1][i-1] = prob
      sample_p_arr[i-1][j-1] = sample_prob
      sample_p_arr[j-1][i-1] = sample_prob

  p_arr = np.array(p_arr)
  sample_p_arr = np.array(sample_p_arr)

  print('inside-outside p:')
  print(p_arr)
  print()

  print('sample p:')
  print(sample_p_arr)
  print()

test('AGACU')
test('AC')
test('GUAC')
test('UUCAGGA')