import math
R = 8.31446

def calculate_mfe(rna, secondary, T):
  energy = {
    'GC': -3,
    'CG': -3,
    'AU': -2,
    'UA': -2,
    'GU': -1,
    'UG': -1
  }

  stack = []
  mfe = 1

  for i, symbol in enumerate(secondary):
    if symbol == '(':
      stack.append(i)
    elif symbol == ')':
      if stack:
        j = stack.pop()
        pair = rna[j] + rna[i]
        mfe *= math.exp(-1*energy.get(pair, 0)/(R*T))
  
  return mfe