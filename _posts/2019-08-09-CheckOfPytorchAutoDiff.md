---
title: Check of automatic differentiation
excerpt: code that verifies that pytorch automatically computed gradient from KL divergence agrees with the formula in the tSNE paper.
type: posts
layout: single
tags: python machine-learning
---

## Comparison of torch automatic gradient with gradient in original python implementation

This code computes the gradient of the tsne KL loss with respect to the positions of the points in the low dimensional space in two different ways: using the formula in the [paper by van der Maaten and Hinton](https://lvdmaaten.github.io/publications/papers/JMLR_2008.pdf) and using the automatic differentiation in pytorch, and verifies that the answers are the same.

The computation using the paper is a (slighly modified) verion of the implementation on [van der Maaten's web page](https://lvdmaaten.github.io/tsne/code/tsne_python.zip).


```python
import torch

# Set up some random data
device='cpu'
P = np.abs(torch.randn(10,10))
P[range(10),range(10)]=0
P = (P + P.t())/2/P.sum()
Y = torch.randn(10,2,requires_grad=True)
Z = Y.clone().detach()
n=P.shape[0] 
no_dims = 2
L=[[(i!=j) for i in range(n)] for j in range(n)]
mask = torch.ByteTensor(L,device=device)
dY=torch.zeros(n,2)

def dist_matrix(Y):
    sum_Y = torch.sum(torch.mul(Y,Y), dim=1)
    num = -2. * (torch.mm(Y,Y.t()))
    num2 = (torch.add(torch.add(num, sum_Y).t(),sum_Y))
    return num2

# Compute the KL loss 
def KL_loss(P,Y,device='cpu'):  
    D = dist_matrix(Y)   
    num2 = 1. / (1. + D)
    numU = torch.masked_select(num2,mask) 
    Q = numU / numU.sum()
    PU = torch.masked_select(P, mask)

    E = (PU*(torch.log(PU/Q))).sum()  
    return E

# compute the gradient using the original code
# note the 4, added to make the agreement exact
def orig_grad(P,Y):
    D = dist_matrix(Y)   
    num2 = 1. / (1. + D)

    num2[range(n),range(n)]=0
    Q = num2 / num2.sum()
    Pa = torch.max(P, torch.tensor([1e-12],device=device))
    Q = torch.max(Q, torch.tensor([1e-12],device=device))
    PQ = P - Q
    for i in range(n):
        M = (PQ[:,i] * num2[:,i]).repeat(no_dims, 1).t()
        B = Y[i,:] - Y
        dY[i, :] = 4*torch.sum(M * B, 0)
    return dY


E = KL_loss(P,Y)
E.backward()

# Doing E.backward() should compute Y.grad, which we save in S
S=Y.grad.clone().detach()
Y.grad.zero_()

# Z is another copy of the Y data, and we compute the gradient by the original method
R = orig_grad(P,Z)
T=S/R
```

The test shows that we are on target.


```python
torch.all(T.isclose(torch.ones(*T.shape))).item()==1
```




    True




```python

```
