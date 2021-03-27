---
title: '替换google字体，加快网站访问速度'
date: 2015-01-07 01:00:00+0800
categories:
    - 博客搭建
tag:
    - google
---

![font log](http://image.fungenomics.com/st.post.2015-01-07-font-face.png)

一段时间以来，这个博客的打开速度慢得出奇！本来只有很少的东西，ping的速度也在120ms左右，算是可以的了，怎么会这样！一开始我没有搞明白问题的根源。今天才突然醒悟到一定是google字体加载的问题！当时这个博客的主题是从[yihui](http://yihui.name/)和[Carl Chen](http://webfrogs.me/)那里抄来的（不用试了，这俩的页面是极难打开的了），因为懒，因为不认识网页语言，尤其是`css`，所以所有`.css`后缀什么的，我一概不去看！所以那几个css就成了我的一个暗区，现在问题很明显了，要修改的地方一定在那！目标锁定之后，那么着手处理吧！

我页面的主要配置都在`style.css`和`home.css`，打开这两文件，搜索`google`，果然发现“fonts.googleapis.com”！那么剩下的就是要先找字体替换方案了，搜索一番后，发现原来360已经把整个的google字库下载下来放在了自己的服务器上了，这样就太好了！改起来也极其简单，只需要直接把所有的`fonts.googleapis.com`替换为`fonts.useso.com`,而且还不用做任何其他的修改，字体也还是原来的字体！

完成后，发现这下页面打开的速度和以前相比根本就是两个网站！不禁感叹：“靠！怎么这么简单！”

*参考 :*   

> <http://www.cgtt.net/1323.html>

欢迎关注我的个人公众号：**helixminer（碱基矿工）**

![helixminer-QRCode](https://static.fungenomics.com/images/2021/03/helixminer-mid-red.png)
