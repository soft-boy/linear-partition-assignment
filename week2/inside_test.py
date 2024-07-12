import unittest
from inside_outside_2 import inside_2

T=1000000000000 # big T makes distribution uniform to mimic total()

class TestInside(unittest.TestCase):
  def test_1(self):
    rna, total = "ACAGU", 6
    
    S_in, _ = inside_2(rna, b=100, T=T)
    self.assertEqual(round(S_in[-1][1]), total)

  def test_2(self):
    rna, total = "AC", 1
    
    S_in, _ = inside_2(rna, b=100, T=T)
    self.assertEqual(round(S_in[-1][1]), total)
  
  def test_3(self):
    rna, total = "GUAC", 5
    
    S_in, _ = inside_2(rna, b=100, T=T)
    self.assertEqual(round(S_in[-1][1]), total)
  
  def test_4(self):
    rna, total = "GCACG", 6
    
    S_in, _ = inside_2(rna, b=100, T=T)
    self.assertEqual(round(S_in[-1][1]), total)

  def test_5(self):
    rna, total = "CCGG", 6
    
    S_in, _ = inside_2(rna, b=100, T=T)
    self.assertEqual(round(S_in[-1][1]), total)
  
  def test_6(self):
    rna, total = "CCCGGG", 20
    
    S_in, _ = inside_2(rna, b=100, T=T)
    self.assertEqual(round(S_in[-1][1]), total)
  
  def test_7(self):
    rna, total = "UUCAGGA", 24
    
    S_in, _ = inside_2(rna, b=100, T=T)
    self.assertEqual(round(S_in[-1][1]), total)
  
  def test_8(self):
    rna, total = "UUGGACUUG", 129
    
    S_in, _ = inside_2(rna, b=100, T=T)
    self.assertEqual(round(S_in[-1][1]), total)

  def test_9(self):
    rna, total = "UUUGGCACUA", 179
    
    S_in, _ = inside_2(rna, b=100, T=T)
    self.assertEqual(round(S_in[-1][1]), total)
  
  def test_10(self):
    rna, total = "GAUGCCGUGUAGUCCAAAGACUUC", 2977987
    
    S_in, _ = inside_2(rna, b=100, T=T)
    self.assertEqual(round(S_in[-1][1]), total)
  
  def test_11(self):
    rna, total = "AGGCAUCAAACCCUGCAUGGGAGCG", 560580
    
    S_in, _ = inside_2(rna, b=100, T=T)
    self.assertEqual(round(S_in[-1][1]), total)

if __name__ == '__main__':
  unittest.main()