from linear_partition import *
import numpy as np
np.set_printoptions(precision=2)

T = 0.01

def test(rna):
  print('p:', rna)

  n = len(rna)
  Q = linear_partition(rna, T=T)
  p, _ = base_pairing_probs(rna, Q, T=T)

  arr = np.zeros((n, n))
  for i in range(1, n+1):
    for j in range(1, n+1):
      prob = p[(i, j)]
      if prob < 0.001:
        prob = 0.0
      arr[i-1][j-1] = prob

  arr = np.array(arr)
  print(arr)
  print()

test('ACAGU')
test('AC')
test('GUAC')
test('UUCAGGA')