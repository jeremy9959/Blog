---
type: posts
layout: single
excerpt: How to combine matplotlib plots into a pdf document 
tags: bioinformatics python
title: PdfPages builds a pdf document from matplotlib figures
---

## How to combine matplotlib figures into a pdf document


This is a useful construct when you want to combine several plots into a single document.

```python
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

x=np.linspace(-3,3,500)

with PdfPages('report.pdf') as pdf:
    for i in range(5):
        fig, ax = plt.subplots(1)
        ax.set_xlim([-3,3])
        ax.plot(x,x**2 + i*x + 1)
        ax.set_title('Plot of x^2+'+str(i)+'*x+1')
        pdf.savefig(fig)
```

[Here is the resulting pdf document](/assets/docs/report.pdf).
