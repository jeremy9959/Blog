---
title: Fixing links with jekyll and gh-pages
excerpt: How to clean up image links for a project page
layout: single
---

I recently moved my blog to a GitHub project page, and this caused trouble with the links
to images. When previewing the blog locally, the jekyll processor wanted to see absolute links
like 

```md
![png](/assets/images/pic.png)
```

but when uploaded to GitHub these links didn't work.

The [solution](https://jekyllrb.com/docs/github-pages/) is buried in  the jekyll documentation
under the heading **Project Page URL Structure**.  The recommendation is:

```
{% raw %}
<!-- For styles with static names... -->
<link href="{{ "/assets/css/style.css" | relative_url }}" rel="stylesheet">
<!-- For documents/pages whose URLs can change... -->
[{{ page.title }}]("{{ page.url | relative_url }}")
{% endraw %}
```



The annoying thing was to go back through all my previous posts and change the links. The following little
python filter/regexp seems to do the trick.

```python
import re
import sys
    
for x in sys.stdin:
	{% raw %}
	m=re.sub(r'(!\[[A-Za-z0-9\-\_ ]+\])\((\S+)\)','<img src=\"{{ \"\\2\" | relative_url }}\">',x.strip())
	{% endraw %}
	print(m)
```

And one another annoyance: to put code with curly braces into a markdown document that later gets processed with jekyll, 
you need to enclose it in "raw/endraw" tags:

{% assign openTag = '{%' %}
```liquid {% raw %}
{% raw %}
{% endraw %}{{ openTag }} endraw %}
```

And as a final remark, I am deeply grateful to [Slaks](https://blog.slaks.net/2013-06-10/jekyll-endraw-in-code/) for showing
how to solve the incredibly annoying problem of getting the endraw tag properly displayed!



