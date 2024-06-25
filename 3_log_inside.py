import unittest
import math
from linear_partition import linear_partition


T=float('inf') # big T makes distribution uniform to mimic total()
phi = lambda x: math.log(x) # log space
logsum = lambda a, b: math.log(math.exp(a)+math.exp(b))
combine = lambda x, y: x+y
log_semiring=(phi, logsum, combine, float('-Inf'), 0)

class TestLogInside(unittest.TestCase):
  def test_1(self):
    rna, total = "ACAGU", 6
    
    Q = linear_partition(rna, b=100, semiring=log_semiring, T=T)
    self.assertEqual(round(math.exp(Q[-1][1])), total)

  def test_2(self):
    rna, total = "AC", 1
    
    Q = linear_partition(rna, b=100, semiring=log_semiring, T=T)
    self.assertEqual(round(math.exp(Q[-1][1])), total)
  
  def test_3(self):
    rna, total = "GUAC", 5
    
    Q = linear_partition(rna, b=100, semiring=log_semiring, T=T)
    self.assertEqual(round(math.exp(Q[-1][1])), total)
  
  def test_4(self):
    rna, total = "GCACG", 6
    
    Q = linear_partition(rna, b=100, semiring=log_semiring, T=T)
    self.assertEqual(round(math.exp(Q[-1][1])), total)

  def test_5(self):
    rna, total = "CCGG", 6
    
    Q = linear_partition(rna, b=100, semiring=log_semiring, T=T)
    self.assertEqual(round(math.exp(Q[-1][1])), total)
  
  def test_6(self):
    rna, total = "CCCGGG", 20
    
    Q = linear_partition(rna, b=100, semiring=log_semiring, T=T)
    self.assertEqual(round(math.exp(Q[-1][1])), total)
  
  def test_7(self):
    rna, total = "UUCAGGA", 24
    
    Q = linear_partition(rna, b=100, semiring=log_semiring, T=T)
    self.assertEqual(round(math.exp(Q[-1][1])), total)
  
  def test_8(self):
    rna, total = "UUGGACUUG", 129
    
    Q = linear_partition(rna, b=100, semiring=log_semiring, T=T)
    self.assertEqual(round(math.exp(Q[-1][1])), total)

  def test_9(self):
    rna, total = "UUUGGCACUA", 179
    
    Q = linear_partition(rna, b=100, semiring=log_semiring, T=T)
    self.assertEqual(round(math.exp(Q[-1][1])), total)
  
  def test_10(self):
    rna, total = "GAUGCCGUGUAGUCCAAAGACUUC", 2977987
    
    Q = linear_partition(rna, b=1000, semiring=log_semiring, T=T)
    self.assertEqual(round(math.exp(Q[-1][1])), total)
  
  def test_11(self):
    rna, total = "AGGCAUCAAACCCUGCAUGGGAGCG", 560580
    
    Q = linear_partition(rna, b=1000, semiring=log_semiring, T=T)
    self.assertEqual(round(math.exp(Q[-1][1])), total)

if __name__ == '__main__':
  unittest.main()