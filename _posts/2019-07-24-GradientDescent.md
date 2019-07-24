---
title: Gradient Descent
layout: single
excerpt: a very simple gradient descent with momentum
---

These notes were motivated by the beautiful article [Why momentum really works](https://distill.pub/2017/momentum/) and
is just a couple of examples of the algorithms for gradient descent and gd with momentum.

We'll look at the case of minimizing a quadratic function $H(w)=\frac{1}{2}w^{t}Aw-b^{t}A$ where
A is a symmetric invertible real matrix and b is a vector.  

We construct $A$ and $b$ as random examples.


```python
Y=np.random.randint(-10,10,size=9).reshape(3,3)
A=np.dot(Y.T,Y)
b=np.random.randint(-10,10,size=3).reshape(3,1)
```


```python
def H(x,A,b):
    """
        evaluate the quadratic function H given by A, b on the vector x.
    """
    return (0.5*np.dot(np.dot(x.T,A),x) - np.dot(b.T,x)).flatten()[0]

def D(x,A,b):
    """        
        evaluate the gradient on the vector x
    """
    return (np.dot(A,x)-b)
```


```python
def descent(alpha, mu, iters, A, b):
    """
        gradient descent with momentum
        alpha: step size
        mu: momentum factor
        iters: number of steps
        A, b: function
    """
    grad=0
    n=b.shape[0]
    x=np.zeros(n).reshape(n,1)
    for i in range(iters):
        grad_new = D(x,A,b)
        grad = mu*grad - alpha*grad_new
        x = x+grad   
    return x
```
