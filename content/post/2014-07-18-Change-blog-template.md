---
title: '博客大改版：添加评论，二维码生成，数学公式的显示，添加分析代码等'
date: 2014-07-18 01:00:00+0800
image: https://image.fungenomics.com/st.post.2014-07-09-github_logo.png
categories:
    - Blog
tags:
    - 网站插件
---

### 目录
----
&nbsp;&nbsp;&nbsp;&nbsp;[**概述**](#introduce)     
&nbsp;&nbsp;&nbsp;&nbsp;[**评论系统由Disqus改成了多说**](#1)    
&nbsp;&nbsp;&nbsp;&nbsp;[**添加google analytics**](#2)     
&nbsp;&nbsp;&nbsp;&nbsp;[**如何搞定数学公式显示问题**](#3)       
&nbsp;&nbsp;&nbsp;&nbsp;[**如何让每个页面自动生成二维码**](#4)   
&nbsp;&nbsp;&nbsp;&nbsp;[**如何添加“返回顶部”按钮“**](#5)   
&nbsp;&nbsp;&nbsp;&nbsp;[**其他**](#6)   
&nbsp;&nbsp;&nbsp;&nbsp;[**最后**](#7)   

<a id='introduce' name='introduce'> </a>

### 概述
----
前几天把博客改版了，在[Jekyll wiki](https://github.com/jekyll/jekyll/wiki/Sites)上爬主题来回爬了好几遍，先是爬了[Linghua Zhang](https://lhzhang.com/),而后又直接爬了[yihui](https://yihui.name/)和[Carl Chen](https://webfrogs.me/)，不过做我这个版面主题的原始作者是[yihui](https://yihui.name/)，随便一提这位是牛人，[统计之都](https://cos.name/)和[COS](https://cos.name/cn/)论坛，以及中国R会议（第一届开始）都是他弄的，年纪轻轻且最近已经出了两本和R相关的书在亚马逊上卖了，恐怕国内（+很多国外）用过R的基本都知道！感兴趣的可以去他的[主页](https://yihui.name/)扒扒。。。扯远了，说回我自己，虽然我也想改得更加不同一些，但一方面暂时还没啥时间，其次，此前未碰过任何与网页制作相关之事，HTML勉强能看，css和js就停留在听说过这两词的地步；然，最重要的是，我也看上了这个版面！所以门面修改一事先缓一缓吧，先在这里谢过各位作者大人！但是话也说回来，以上只要有任意一位作者介意，我也只能作罢，重新扒过。

OK，既然现在外在的部分还无力去动，那接下来，我说一下自己所做的一些内在改变。

<a id='1' name='1'> </a>

### 1. 评论系统由Disqus改成了多说
----

![多说](https://www.lagou.com/upload/webproduct/ff80808142c5ed7f0142c6bb08cc12fe.png)

我倒不是排斥Disqus，一开始我用的就是它，Disqus，国际化，版面简洁，管理容易，主页也生动，cool，我很喜欢！一个账号说遍天下，当然这一点上多说也一样。换掉它根本的原因还是在于伟大的‘墙’，Disqus只具有分享到Facebook和Twitter的功能，而这两货正常途径咱是上不了的。也罢，在国内的话多说用起来确实要更友好些，所以这次改版就重新选择了和国内社交网络联系在一起的多说。

这里我说一下自己是如何添加的。

（1）如果还没有绑定多说，那第一步需要做的就是到[多说](https://duoshuo.com/)上绑定。多说是不需要做任何注册的，只要有QQ，微博，百度等账号就行，这一点相当方便，不然我又得增加一个账号，我都已经记不得自己在网络上到底注册过多少个账号了，恐怕各型各色几百个都有了！ 

（2）尔后，登陆多说主页点击“我要安装”，按提示走就行了。在`站点地址`这一栏填入自己的域名，再填上其他信息，按"创建"

（3）创建了之后，在侧边框找到选择`工具`，选择`通用代码`，会看到多说提供了一段相关的代码，我们要做的事情很简单，直接把它复制并贴到你自己网页的代码中，然后，就没有然后了。。。诶诶，慢着呀！那应该贴进哪一个的网页里头啊？我这简单说明一下，比如，`_layout`（生成网站页面的相关网页绝大数都放这，jekyll会自己来加载）目录下我们一般都会创建`defualt.html`, `page.html`, `post.html`等这么几个文件。我假设`post.html`就是你用于博客的，而你也只想博客的页面才加多说评论框，那么代码就应该贴进`post.html`。虽然你可以放在<body>与</body>之间的任何位置，但考虑到页面主要内容的加载速度，最好还是把它靠后放，评论模块要那么急着加载干嘛，最好就只是在</body>之前，下面我给了个示意。

添加多说的代码结构。

```html
<html>
<head>
<title></title>
</head>
<body>
...
<!—多说代码-->
{ % include duoshuo.md % } <!-- 实际使用中 “{” 和 “%”，以及“%”和”}“ 之间的空格需要去掉-->
</body>
</html>
```

如果除了post.html，你还打算将它添加到其他的页面中，也没问题，操作和上面类似。但整块整块的代码复制来复制去的，总觉得很麻烦，一旦想修改还得一个个去改，太吃力了！！我自己更喜欢的一种方式就是直接把多说的代码写到一个单独的文件中，比如就叫`duoshuo.md`，并放在`_include`目录下（这个目录按照jekyll的建议，就是让我们用来存放模块插件的），然后在需要多说的页面中直接调用`{ % include duoshuo.md % }`就行(本文为了显示的原因在`{`和`%`以及`%`和`}`之间加了一个空格，正式使用的时候注意去掉，下同)，这就犹如一个函数一般，修改也只需要改duoshuo.md这个文件，相当方便。

在设置`<div class="ds-thread" data-thread-key="请将此处替换成文章在你的站点中的ID" data-title="请替换成文章的标题" data-url="请替换成文章的网址"></div>`，这几个值的时候，我参考了这篇[博文](https://www.leejianyang.com/2014/05/25/duoshuo_tutorial/)，它讲的很清楚，不过里面写成data-thread-key=”<%= page.path %>”，data-title=”<%= page.title %>”和data-url=”<%= page.permalink %>”的形式在我这并不能成功，不知道是否与我用markdown格式编辑有关。所以后来我改成`<div class="ds-thread" data-thread-key="{ { page.path } }" data-title="{ { page.title } }" data-url="{ { page.permalink } }"></div>`就成功了（本文为了显示的原因在`{`和`{`以及`}`和`}`之间加了一个空格，正式使用的时候注意去掉，下同）。

<a id='2' name='2'> </a>  

### 2. 添加google analytics
----

![GA](https://www.proyectosbds.com/blog/wp-content/uploads/2012/11/google-analytics-impementation.jpg)
对于自己的网站，一个关心的事情就是网站的浏览量如何？是通过什么途径以什么方式流入流出的？别的先不多说，光是能时刻知道自己的网站在哪里、被用什么操作系统、多少人浏览、网页之间的流入流出是怎么样的等状态信息的本身就是一个很cool的事情，这跟玩游戏一样，相信没多少人会愿意把自己角色的血条和蓝条隐藏掉，然后在那瞎玩。所以我就开始琢磨着应该怎么搞定这样一个事情。后来知道了可以用[google-analytics](https://www.google.com/analytics/)帮助记录和分析，同样需要登陆google账号，填入网站域名，获得对应该域名的google analytics网络跟踪代码，然后又是复制粘贴。改过一次评论系统之后，我就发现，这东西都是大同小异的，就那么几下板斧。当然了我还是把它单独写成一个模块放在`_include`下供调用。代码如下：

```html
<!-- Google Analysis For website tracking -->
<script>
  (function(i,s,o,g,r,a,m){i['GoogleAnalyticsObject']=r;i[r]=i[r]||function(){
  (i[r].q=i[r].q||[]).push(arguments)},i[r].l=1*new Date();a=s.createElement(o),
  m=s.getElementsByTagName(o)[0];a.async=1;a.src=g;m.parentNode.insertBefore(a,m)
  })(window,document,'script','//www.google-analytics.com/analytics.js','ga');

  ga('create', 'UA-52659904-2', 'auto');
  ga('require', 'displayfeatures'); 
  ga('send', 'pageview');
</script> 
<!-- Google Analysis end -->
```

这一次我希望全站跟踪，`default.html`是我所有页面都会添加的网页，所以这个代码就放在`default.html`中了。关于GA代码应该放在文件中的哪个位置比较适合，我还是做了一下考虑的，参考了[这篇文章](https://bluewhale.cc/2010-07-19/google-analytics-add-tracking-code.html),按照异步跟踪的方式添加：

```html
<html>
<head>
<title></title>
...
<!—GA异步追踪代码-->
{ % include googleAnalysis.md % } <!-- 实际使用中 “{” 和 “%”，以及“%”和”}“ 之间的空格需要去掉-->
...
</head>
<body>
</body>
</html>
```


__注意:__ *在国内使用google analytics 并不是一个明智的做法，它已被封，这会大大降低你的网页访问速度，体验极差，建议还是使用百度分析.*

<a id='3' name='3'> </a>

### 3. 如何搞定数学公式显示问题
----

![mathJax](https://media.journals.elsevier.com/content/images/main/elsevier-continues-as-mathjax-12092220.jpg)

makrdown是个好东西，jekyll+markdown做网站很容易！但麻烦的是，markdown不能支持LaTeX，我写个数学公式在那上面，它解析不了，只是原封不动的显示在那。。这一点真让人捉急！咋办？又是一番的google，后来还是找到了个好办法——[MathJax](https://www.mathjax.org/)，它是一个数学公式显示引擎，能够把LaTeX编辑的公式在网页上显示出来，不过它是在线解析的，公式如果比较多的话，显示速度会稍慢。要用MathJax，需要先把`_config.yml`中的makrdown解析器换成`kramdown`，不然用不了。我同样也是写成模块`mathJax.md`，然后由post.html调用(估计我也不会在除了博客之外的地方用到它)，为了加载速度，也放的靠后一些，刚好在多说模块上面，下面是我用的mathJax.md模块的代码，当然了，若有需要你也可以使用：

```html
<script type="text/x-mathjax-config">
  MathJax.Hub.Config({ tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]} });
</script>
<script type="text/javascript"
  src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML">
</script>
```

成功添加了之后只需在公式前后用\$\$围起来，就可方便的编写公式了。如：

```bash
$$
\rho_i=\sum_j\chi(d_{ij}-d_c)
$$
```

显示的结果就是：
$$
\rho_i=\sum_j\chi(d_{ij}-d_c)
$$

<a id='4' name='4'> </a>

### 4. 如何让每个页面自动生成二维码
----
![QR](https://image.fungenomics.com/st.post.2014-07-18-Figure4-QRcode.jpg)

二维码的作用其实也不必做多解释，最重要的就是方便。我之前并不知道可以用Google API为网站自动生成二维码一事，直到我成功了之后才恍然一悟，真是隔行如隔山。一开始我都是暴力解决：找个能生成二维码的网站，把自己每一篇博客的网址贴进去，点击生成二维码，然后再把这个生成的二维码图片下载下来，接着将它上传至图床，最后再将图片的链接添加到每篇博客的最后！这实在是。。。这个过程就算只是这样说起来都觉得十分费事，可以想象操作起来该有多费劲，而且一旦网站发生调整，就得重来一遍，实在低效！所以我一直寻思着该如何做才能让网页自己去产生二维码，一定要将它自动化！后来注意到了这个[QR生成器](https://www.the-qrcode-generator.com/)，并且受到它`share`按钮中内容的启发，我当时就在想为什么它这样的一段代码（如下）在插入到网页中之后就能出来一个二维码呢？（当然这个二维码是指向它自己网站的）

```html
<!--QR生成器的share代码-->
<a href="https://www.the-qrcode-generator.com/"><img src="https://chart.apis.google.com/chart?chs=200x200&amp;cht=qr&amp;chld=|1&amp;chl=" alt="QR Code" /></a>
```

在仔细观察后发现代码中url地址中后面的`chs=200x200&amp;cht=qr&amp;chld=|1&amp;chl=`看起来很像某些参数，貌似是可以自己设置的，若真如此的话，又该怎么做呢？这次得感谢[DIVCSS5](https://www.divcss5.com/template/m504.shtml)这篇博文中的那一段代码结构，它完全让我明白了自己该如何改装上面的QR生成器代码让它能为我所用，给每一个网页都自动根据自己的网址产生正确的二维码，改了的代码如下：

```html
<!-- 网页自动生成二维码的代码 -->
<script type="text/javascript">  
  thisURL  = document.URL;  
  strwrite = "<p align='center'><img src='https://chart.apis.google.com/chart?chs=120x120&amp;cht=qr&amp;chld=|1&amp;chl=" + thisURL + "' width='120' height='120' alt='QR Code'/>（传送门）</p>";
  document.write( strwrite ); 
</script>
```

这一段代码本身是独立的，不依赖于任何特定的网站拥有者，而且光看代码我们也能想象得到，这里头调用的就是Google API，所以如果有需要，也欢迎大家直接复制粘贴到自己的页面代码中。当然了，我还是把它独立写在了一个文件——`qrCodeGenerate.md`中并放在了`_include`目录下，然后再调用，就跟上面添加多说和google跟踪代码一个思路，不过这个倒没什么位置限制，就看各自喜欢放哪就放哪，只要是在<body>与</body>之间就行，但我还是建议放在博文内容加载模块之后，这样可以确保博文内容先加载。

__注意:__ *后来在国内发现上面这个基于google的二维码，不但没用还大大影响了页面的加载速度！原因大家也都知道，具体的修改请参看这篇[修改二维码生成工具](/2014/09/19/Change-QRCode-Generator.html)博文。*

<a id='5' name='5'> </a>

### 5. 如何添加“返回顶部”按钮
----

![Backtop](https://85ryan.com/wp-content/uploads/2013/09/simple-jquery-back-to-top-button.jpg)

"返回顶部"按钮（本博客右下角）的添加完全参照[这个教程](https://liberize.me/tech/jekyll-add-back-to-top-button.html)，作者写的相当详细！按照文章的信息我把'二'中的代码写入backtop.js，'三'中的信息写入backtop.css。然后将它们放在defualt.html中的<head>与</head>之间调用，示意如下：

```html
<html>
<head>
...
<link rel="stylesheet" href="/media/css/backtop.css">  <!-- Back Top -->
<script type="text/javascript" src="/media/js/backtop.js"></script>  <!-- Back Top -->
...
</head>
<body>
</body>
</html>
```

<a id='6' name='6'> </a>

### 6. 其他
----

（1）关于主页上的“[FORUM](/forum)”，本来有打算做成类似于论坛形式可以让大家就某个话题来进行讨论，但是发现要实现这个功能难度很大工作量不小！限于我自己能力有限，这个也得慢慢来，暂时将就变成了贴“鸡汤文”的板子了，o(╯□╰)o！

（2）google网站站长工具，上传了一份google提供的验证网页（它仅是验证之用），成功之后虽然可以将它删掉，但google官网的意思是建议留下。

（3）图床的选择。图床对于独立博客来说是一个重要又头疼的事情，目前我直接用[点点博客](https://www.diandian.com/)的图片功能，使用外链。它的好处是免费，稳定，访问速度也快，要是真有一天点点不允许这样用的话，那我也不怕，实在不行就用七牛(__注: 已改用七牛__ )，或者想办法用Github作为图床。

<a id='7' name='7'> </a>

### 7. 最后
----
OK！唠叨了不少总算把这篇文章写完了。接下来与博客搭建相关的内容先暂放一边了，我要回归主业的更新了。


欢迎关注我的个人公众号：**helixminer（碱基矿工）**

![helixminer-QRCode](https://static.fungenomics.com/images/2021/03/helixminer-mid-red.png)
