# 1-nonterminal Nussinov
```
S -> \epsilon                    p = 1

# S -> S .
S -> S c | S g | S a | S u       p = 1

# S -> S(S)
S -> S c S g | S g S c           p = e ^ -3/RT
S -> S a S u | S u S a           p = e ^ -2/RT
S -> S g S u | S u S g           p = e ^ -1/RT
```
