---
title: '为什么uBAM迟迟无法流行起来'
date: 2018-06-09 01:00:00+0800
image: https://static.fungenomics.com/images/2021/03/bam-ex-20210327230528877.jpg
categories:
    - 生物信息
    - 基因组学
tags:
    - BAM


---



**uBAM就是非比对的BAM文件，**fastq可以通过picard这个工具将其转为这个格式。

它有不少优于fastq格式的地方，比如：同一个read的数据都在同一行；拓展性强，可以添加丰富的metadata；方便维护，同一个样本的测序数据甚至可以只通过一份uBAM来存储就行了等。

从我第一次知道uBAM的时候算起，已经过去4+年了。我也曾非常看好它，认为它必是以后存储下机数据的标准。**然而时隔多年，很奇怪，uBAM这么好（GATK也都一直支持着这个格式），为什么还是迟迟不见它流行起来呢？**

目前，使用uBAM格式的单位据我所知也仅仅只是一些比较大型的研究机构，比如美国的Broad Institute和英国的Sanger会采用它来存储下机数据。

这段时间思考下来觉得可能有以下几个原因，与诸位共享：

1. **BAM是“笨重”的**，它并不是文本文件，你无法直接通过文本工具打开它查看具体内容。只能通过第三方工具或者专门的SAM/BAM程序包（或者API）来实现对它的操作。这对许多不熟悉这一处理方式的研究者来说，会带来很多麻烦。这等于是直接提高了操作这个文件的门槛，从这一点看使用体验确实远不如fastq；

2. **主流工具还不完全支持**，除了samtools和与它相关的少量工具，并没有太多其他的工具直接支持在命令行操作BAM；

3. **BAM文件的空间占比并不比压缩了的fastq小很多**，优势有限；

4. 底层IO效率方面，实际上也是文本格式的fastq（或者gzip压缩的fastq）要高于BAM。

从uBAM的这个现象，或许也侧面折射出了一些关于产品设计（或者方案设计）的问题。关于这个问题，我看到了三个地方，欢迎大家拍砖：

第一、体验。一个产品或者方案要流行起来，除了解决需求之外，对 **使用体验的关注度要高于技术的先进性和产品本身的完备性；**

第二、先发优势。时间一旦落后了（比如fastq早于uBAM很多年），用户习惯的更改需要有完备的技术解决工具来支持，降低切换成本，甚至实现无痛切换，从而最大程度的保留新产品的优势；

第三、看似简单的事物越是难以被取缔。fastq格式是一个存储测序数据极为简单、简明的数据格式，它只包含所有必须包含的内容，而且目标明确，就是序列ID、测序数据和质量值，它们都是必不可缺的信息，再多无用，似乎已是极致。

* * *

推荐阅读

* [从零开始完整学习全基因组测序数据分析：第2节 FASTA和FASTQ](http://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798261&idx=1&sn=48d277f96ac65ed66f2e5d06f11b5f14&chksm=83c1d469b4b65d7fb634059356252989b03a3e3e34e6dfe12b88ada6d8a52ab9b569f088c839&scene=21#wechat_redirect)

* [GATK4.0和全基因组数据分析实践（上）](http://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798425&idx=1&sn=ae355ed362848578e5c853413f23dfd7&chksm=83c1d505b4b65c13124c9acd210356c4364ec9f5498bbd16fa4475be29811213abb64ea9720f&scene=21#wechat_redirect)

* [GATK4.0和全基因组数据分析实践（下）](http://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798455&idx=1&sn=67a7407980a57ce138948eb46992b603&chksm=83c1d52bb4b65c3dde31df94e9686654bf616166c7311b531213ebf0010f67a32ce827e677b1&scene=21#wechat_redirect)

* [我是解螺旋的矿工，我热爱生命科学](http://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798476&idx=1&sn=52cd7f2d86d77410d538ab9edbec8df7&chksm=83c1d550b4b65c468425dfc367c23c4da80ebfe7060df2bd314354a3b4690c9425640153ba51&scene=21#wechat_redirect)

* [该如何自学入门生物信息学](http://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798366&idx=1&sn=b545fcea7f82839fa87e9d9e472d1e72&chksm=83c1d4c2b4b65dd4843250c307969ada96c4039f4f528c034620d25b78d8beba2f9cf924bb8a&scene=21#wechat_redirect)

***

欢迎关注我的个人公众号：**helixminer（碱基矿工）**

![helixminer-QRCode](https://static.fungenomics.com/images/2021/03/helixminer-mid-red-20210327230505192-20210327230529031.png)

***

这是知识星球：『解螺旋技术交流圈』，是一个我与读者朋友们的私人朋友圈。我有9年前沿而完整的生物信息学、NGS领域的工作经历，在该领域发有多篇Nature级别的科学文章，我也希望借助这个知识星球把自己的一些微薄经验分享给更多对组学感兴趣的伙伴们。

自从星球正式运行以来，已经过去了6个月，星球的成员也已经超过220人了。所分享的主题超过了500个，回答的问题超过了140个，精华70个。我在知识星球上留下的文字估计也已经超过10万字，加上大家的就更多了，相信接下来星球的内容一定还会不断丰富。另外，上周获得了知识星球官方评选的“最优质星球”优秀奖。

这是知识星球上 **第一个真正与基因组学和生物信息学强相关的圈子。**我希望能够借此营造一个高质量的组学知识圈和人脉圈，通过提问、彼此分享、交流经验、心得等，彼此更好地学习生信知识，提升基因组数据分析和解读的能力。

在这里你可以结识到全国优秀的基因组学和生物信息学专家，同时可以分享你的经验、见解和思考，有问题也可以向我提问和圈里的星友们提问。

知识星球邀请链接：[「解螺旋技术交流圈」](https://wx.zsxq.com/mweb/views/joingroup/join_group.html?group_id=518881585444&secret=vcdvs4rdpst7stq4wcvqmlwvogc0ssbn&user_id=28821152428221)