# inside
```
for j = 1 ... n:
  for i s.t. [i, j-1] in \in(S):
    
    \in(S_i,j) += \in(S_i,j-1) * (e^-delta(j)/RT)

    if i-1 and j pair:
      for k s.t. [k, i-2] in \in(S):
        \in(P_i-1,j) = \in(S_i,j-1) * (e^-xi(i-1,j)/RT)
        \in(S_k,j) += \in(S_k,i-2) * \in(P_i-1,j)
  
  beamprune()
```
# outside
```
for j = n ... 1:
  for i s.t. [i, j-1] in \in(S):
    
    \out(S_i,j-1) += \out(S_i,j) * (e^-delta(j)/RT)

    if i-1 and j pair:
      for k s.t. [k, i-2] in \in(S):

        \out(S_k,i-2) += \out(S_k,j) * \in(P_i-1,j)
        \out(P_i-1,j) += \out(S_k,j) * \in(S_k,i-2)
        \out(S_i,j-1) += \out(P_i+1,j) * (e^-xi(i-1,j)/RT)

      prob_i-i,j = \out(P_i-1,j) * \in(P_i-1,j) / total
  
  beamprune()
```