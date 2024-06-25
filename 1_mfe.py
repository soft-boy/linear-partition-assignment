import unittest
from linear_partition import linear_partition
from calculate_mfe import calculate_mfe

T=1
phi = lambda x: x # phi is the homomorphic function from R+ to the new set
summarize = lambda x, y: max(x, y)
combine = lambda x, y: x*y
semiring=(phi, summarize, combine, 0, 1)

class TestInside(unittest.TestCase):
  def test_1(self):
    rna, secondary = "ACAGU", '((.))'
    
    Q = linear_partition(rna, b=100, semiring=semiring, T=T)
    mfe = calculate_mfe(rna, secondary, T)
    self.assertEqual(round(Q[-1][1]), round(mfe))

  def test_2(self):
    rna, secondary = "AC", '..'
    
    Q = linear_partition(rna, b=100, semiring=semiring, T=T)
    mfe = calculate_mfe(rna, secondary, T)
    self.assertEqual(round(Q[-1][1]), round(mfe))
  
  def test_3(self):
    rna, secondary = "GUAC", '(())'
    
    Q = linear_partition(rna, b=100, semiring=semiring, T=T)
    mfe = calculate_mfe(rna, secondary, T)
    self.assertEqual(round(Q[-1][1]), round(mfe))
  
  def test_4(self):
    rna, secondary = "GCACG", '().()'
    
    Q = linear_partition(rna, b=100, semiring=semiring, T=T)
    mfe = calculate_mfe(rna, secondary, T)
    self.assertEqual(round(Q[-1][1]), round(mfe))

  def test_5(self):
    rna, secondary = "CCGG", '(())'
    
    Q = linear_partition(rna, b=100, semiring=semiring, T=T)
    mfe = calculate_mfe(rna, secondary, T)
    self.assertEqual(round(Q[-1][1]), round(mfe))
  
  def test_6(self):
    rna, secondary = "CCCGGG", '((()))'
    
    Q = linear_partition(rna, b=100, semiring=semiring, T=T)
    mfe = calculate_mfe(rna, secondary, T)
    self.assertEqual(round(Q[-1][1]), round(mfe))
  
  def test_7(self):
    rna, secondary = "UUCAGGA", '(((.)))'
    
    Q = linear_partition(rna, b=100, semiring=semiring, T=T)
    mfe = calculate_mfe(rna, secondary, T)
    self.assertEqual(round(Q[-1][1]), round(mfe))

if __name__ == '__main__':
  unittest.main()
