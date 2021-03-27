---
title: '修改二维码生成工具'
date: 2014-09-19 01:00:00+0800
image: http://image.fungenomics.com/fg.post.qrcode.png
categories:
    - 博客搭建
tags:
    - 网站插件
    - google
---

本来博客在改版后，短时间内，我是不想做修改了的，回国后，没办法，原来的用于生成二维码的代码必须得换！根本就是失效了。原因也是很明朗的，二维码生成代码用的API是google的，稳定，强大......但在墙面前显得跟土一样。有时候我总觉得，是不是因为有了墙，大家的技术反而都要强那么一点了——魔高一尺道高一丈！只是有些时候也不得不拜倒了。总之，为了速度，为了大家的良好体验，换了！google的跟踪代码，我也取消了，换成了百度，不然网页死活打不开的，这是题外话了。

总之货比三家，最后我做了个比较大众化的选择——联图网的服务接口！说是接口，听起来似乎很高大上的样子，其实也就是一小段代码罢了，无需紧张。下面就是主题了。

API接口调用代码如下：

```html
<!-- 网页自动生成二维码的代码 -->
<script type="text/javascript">  
  thisURL  = document.URL;  
  strwrite = "<p align='center'><img src='http://qr.liantu.com/api.php?w=120&m=2&text=" + thisURL + "' alt='QR Code'/>（传送门）</p>"
   document.write( strwrite ); 
</script>
```

还和上次[博客改版](/2014/07/18/Change-blog-template.html)中提到的一样，把上面这段代码完全替换掉原来`qrCodeGenerate.md`中的内容就ok了。

为了方便大家看的明白，这里需要对上面这段代码中的一些重要参数做些说明：

1. `http://qr.liantu.com/api.php?` 是联图网的QR code API接口地址；   
2. w：二维码图片宽度，根据需要自行修改，如`w=120`，表示生成一个120×120像素的正方形二维码图片；     
3. m：二维码静区（外边距），一般设为比较小的数，我这里设置为2；    
4. text：二维码数据信息， 我这里设置为对应页面的URL地址；     
5. alt：二维码图片描述， 这个就没什么讲究了。 

总体的效果请看文章末尾的二维码图片。

*参考资料：*

> [免插件自动生成wordpress文章二维码图片](http://uuxn.com/wordpress-articles-qr)  很感激这篇文章给我的启发。


欢迎关注我的个人公众号：**helixminer（碱基矿工）**

![helixminer-QRCode](https://static.fungenomics.com/images/2021/03/helixminer-mid-red.png)




