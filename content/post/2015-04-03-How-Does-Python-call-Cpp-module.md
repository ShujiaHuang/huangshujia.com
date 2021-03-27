---
title: '在Python中调用C++模块'
date: 2015-04-03 01:00:00+0800
categories:
    - 信息图
    - 编程技术
tags:
    - C++
    - Python
---


在Python中成功实现了对原来C++代码模块的复用！这个好处多多，Python写得快，C++跑得快，那就是既快又快了！方法很简单，以至于我能够用一张截图记录下整个过程（点击图片看大图）！

[![Python-call-C++-Module](http://image.fungenomics.com/st.post.2015-04-03-Python-call-Cpp-Module.png)](http://image.fungenomics.com/st.post.2015-04-03-Python-call-Cpp-Module.png)

其实，注意到，必须在原来的C++代码后面添加extern "C"来辅助（C则不需要，这也是与复用C代码时最大的不同点），不然Python在调用这个构建后的动态链接库时是找不到原来的方法或者函数的，说到底还都是因为当前Python的设定中只能调用C函数，而不能直接调用C++的方法，因此extern "C"封装的函数必须是C风格的，也即就是说，函数中的所有参数（输入参数和返回参数）都不能包含任何非C定义的类型（包括不能包含只属于C++的类型）。


欢迎关注我的个人公众号：**helixminer（碱基矿工）**

![helixminer-QRCode](https://static.fungenomics.com/images/2021/03/helixminer-mid-red.png)



