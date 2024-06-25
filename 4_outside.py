import unittest
from linear_partition import linear_partition, base_pairing_probs

T=1000000000000 # big T makes distribution uniform to mimic total()
phi = lambda x: x # phi is the homomorphic function from R+ to the new set
plus = lambda x, y: x+y
mult = lambda x, y: x*y
semiring=(phi, plus, mult, 0, 1)

class TestOutside(unittest.TestCase):
  def test_1(self):
    rna, total = "ACAGU", 6
    
    Q = linear_partition(rna, b=100, semiring=semiring, T=T)
    _, Q_hat = base_pairing_probs(rna, Q, T=T)
    
    self.assertEqual(round(Q_hat[(1, 0)]), round(Q[-1][1]))
    self.assertEqual(round(Q_hat[(1, 0)]), total)

  def test_2(self):
    rna, total = "AC", 1
    
    Q = linear_partition(rna, b=100, semiring=semiring, T=T)
    _, Q_hat = base_pairing_probs(rna, Q, T=T)
    
    self.assertEqual(round(Q_hat[(1, 0)]), round(Q[-1][1]))
    self.assertEqual(round(Q_hat[(1, 0)]), total)
  
  def test_3(self):
    rna, total = "GUAC", 5
    
    Q = linear_partition(rna, b=100, semiring=semiring, T=T)
    _, Q_hat = base_pairing_probs(rna, Q, T=T)
    
    self.assertEqual(round(Q_hat[(1, 0)]), round(Q[-1][1]))
    self.assertEqual(round(Q_hat[(1, 0)]), total)
  
  def test_4(self):
    rna, total = "GCACG", 6
    
    Q = linear_partition(rna, b=100, semiring=semiring, T=T)
    _, Q_hat = base_pairing_probs(rna, Q, T=T)
    
    self.assertEqual(round(Q_hat[(1, 0)]), round(Q[-1][1]))
    self.assertEqual(round(Q_hat[(1, 0)]), total)

  def test_5(self):
    rna, total = "CCGG", 6
    
    Q = linear_partition(rna, b=100, semiring=semiring, T=T)
    _, Q_hat = base_pairing_probs(rna, Q, T=T)
    
    self.assertEqual(round(Q_hat[(1, 0)]), round(Q[-1][1]))
    self.assertEqual(round(Q_hat[(1, 0)]), total)
  
  def test_6(self):
    rna, total = "CCCGGG", 20
    
    Q = linear_partition(rna, b=100, semiring=semiring, T=T)
    _, Q_hat = base_pairing_probs(rna, Q, T=T)
    
    self.assertEqual(round(Q_hat[(1, 0)]), round(Q[-1][1]))
    self.assertEqual(round(Q_hat[(1, 0)]), total)
  
  def test_7(self):
    rna, total = "UUCAGGA", 24
    
    Q = linear_partition(rna, b=100, semiring=semiring, T=T)
    _, Q_hat = base_pairing_probs(rna, Q, T=T)
    
    self.assertEqual(round(Q_hat[(1, 0)]), round(Q[-1][1]))
    self.assertEqual(round(Q_hat[(1, 0)]), total)
  
  def test_8(self):
    rna, total = "UUGGACUUG", 129
    
    Q = linear_partition(rna, b=100, semiring=semiring, T=T)
    _, Q_hat = base_pairing_probs(rna, Q, T=T)
    
    self.assertEqual(round(Q_hat[(1, 0)]), round(Q[-1][1]))
    self.assertEqual(round(Q_hat[(1, 0)]), total)

  def test_9(self):
    rna, total = "UUUGGCACUA", 179
    
    Q = linear_partition(rna, b=100, semiring=semiring, T=T)
    _, Q_hat = base_pairing_probs(rna, Q, T=T)
    
    self.assertEqual(round(Q_hat[(1, 0)]), round(Q[-1][1]))
    self.assertEqual(round(Q_hat[(1, 0)]), total)
  
  def test_10(self):
    rna, total = "GAUGCCGUGUAGUCCAAAGACUUC", 2977987
    
    Q = linear_partition(rna, b=100, semiring=semiring, T=T)
    _, Q_hat = base_pairing_probs(rna, Q, T=T)
    
    self.assertEqual(round(Q_hat[(1, 0)]), round(Q[-1][1]))
    self.assertEqual(round(Q_hat[(1, 0)]), total)
  
  def test_11(self):
    rna, total = "AGGCAUCAAACCCUGCAUGGGAGCG", 560580
    
    Q = linear_partition(rna, b=100, semiring=semiring, T=T)
    _, Q_hat = base_pairing_probs(rna, Q, T=T)
    
    self.assertEqual(round(Q_hat[(1, 0)]), round(Q[-1][1]))
    self.assertEqual(round(Q_hat[(1, 0)]), total)

if __name__ == '__main__':
  unittest.main()