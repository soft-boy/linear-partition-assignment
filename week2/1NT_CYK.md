# 1-nonterminal CYK
```
S -> c | g | a | u       p = 1
S -> S S                 p = 1

# S -> (S)
S -> c S g | g S c       p = e ^ -3/RT
S -> a S u | u S a       p = e ^ -2/RT
S -> g S u | u S g       p = e ^ -1/RT
```
