---
layout: single
title: Some emacs/beamer notes from today
excerpt: a few notes on beamer formatting (code, colors)
---

## A few notes on beamer from this  morning

### How to include python code in a beamer presentation

With thanks to [felix11h.github.io](http://felix11h.github.io/blog/latex-beamer-minted).

- You need the minted package
- You need to give the shell escape option to TeX
- you need to label the frame environments as fragile
- you need to have pygmentize installed and accessible (so perhaps activate the right venv)


```
\documentclass{beamer}
\usepackage{minted}
\newminted{python}{fontisze=\scriptsize,
%                  linenos,
                   numbersep=8pt,
%                  frame=lines,
                   gobble=4,
                   framesep=3mm}
\begin{document}
\begin{frame}[fragile]
\frametitle{Sample Code}
\begin{pythoncode}                                                                                                     
	import numpy as np
	def f(x):
		return x**2
\end{pythoncode}
\end{frame}
\end{document}

%%% Local Variables:
%%% mode: latex
%%% Tex-master: t
%%% TeX-command-extra-options: "-shell-escape"
%%% End:
```


