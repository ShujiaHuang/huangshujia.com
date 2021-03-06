---
title: '样本量重要，还是测序深度重要? 生物信息工程师可以分为多少种类型? |《解螺旋技术交流圈》精华第3期'
date: 2018-06-17 01:00:00+0800
image: https://static.fungenomics.com/images/2021/03/helixminer-club.png
categories:
    - 生物信息
    - 基因组学
    - 《解螺旋技术交流圈》精华
tags:
    - 样本量
    - WES
    - 职业发展


---



今天，继续把发在“解螺旋技术交流圈”的部分主题整理出来，分享给你。

## 1. 请问对于同一份BAM文件使用samtools depth和用samtools mpileup跑出来的位点的depth有何差异？

你会注意到这个差异，应该是由于你所用的是Pair-End（PE）测序的数据吧，如果是SE数据，差异其实很小。对于PE测序数据主要有两个地方的差异：

![](https://static.fungenomics.com/images/2021/03/samtools-mpileup-20210327230739795.png)

<p align="center"><a>samtools mpileup</a></p>

（1）第一个差异，对于PE数据，mpileup默认会把不正常比对的PE Read（比如read1和read2的比对位置彼此间的距离超过插入片段长度的波动范围或者read1与read2有一条没有比对上）先排除掉再做计算，但samtools depth则不会，depth默认不做任何过滤，只要比上就算。这也是我们会看到samtools depth计算的覆盖深度往往都高于mpileup的最主要原因。如果要让两者一致，可以在mpileup中加上 -A 参数，强制留下不正常的PE比对结果即可；

（2）它们之间的第二个差异是，在默认情况下，mpileup还会过滤掉测序质量值低于13的碱基，depth默认不过滤。

虽然调整一下参数就可以保证两者一样。但我并不建议这么做，虽说mpileup这里得到的是高质量的覆盖深度，但是说到底它和samtools depth的目的还是不同的。

此外，如果要更好地计算比对数据的覆盖深度和覆盖度的话，samtools depth虽然能够胜任，但是功能还是比较单一，而且由于每个位点都会输出，导致结果文件总是很巨大，我还是比较推荐使用bedtools2来完成，如下图，它的功能和输出形式要更加丰富。

![](https://static.fungenomics.com/images/2021/03/bedtools-20210327230739891.png)

<p align="center"><a>bedtools2计算基因组覆盖度的不同模式</a></p>

## 2. 为什么WES的数据无法使用VQSR进行变异质控？

其实不只是WES，还包括很多小panel的数据，如果样本量比较少的话基本都无法使用VQSR进行变异的质控。其原因就在VQSR的原理上。

VQSR的核心原理是利用机器学习算法构造一个区分“好”变异和“坏”变异的分类器。这个分类器在GATK中是通过GMM模型来构造的，它在构造的时候并不是盲目地使用所有数据来进行构造，而是挑出和已知的变异集合Overlap的位点（通常是HapMap数据集）——并分配相应的可信度权重来进行训练。

基于群体遗传的原理，这些已知且被严格验证的变异（如HapMap数据）会被认为是更加靠谱的变异，因此在初始化的时候先把它们当作是“好”的——也就是正确的变异。这个初始变异集很重要，然后利用这些好变异训练一个区分好变异的GMM，接着对全部数据进行打分，再把评分最低的那些拿出来，构成一个最不像正确变异的集合，用来构造一个区分坏变异的GMM，用来专门识别坏变异。最后同时用好和坏的GMM再一次同时对变异进行打分，看每个变异更像谁，就能够评判出这个变异可信的质量值了。越靠近好的GMM，质量就越高，这就是VQSR过滤的大致原理（如下图）。

![ ](https://static.fungenomics.com/images/2021/03/vqsr_model-20210327230739966.png)

<p align="center"><a>VQSR区分好变异和坏变异的分类器</a></p>

为了得到理想好的结果，VQSR在进行模型训练的时候就有一个最低可用位点数目的要求——通常是好和坏变异可供训练的数目必须超过5000个，如果Overlap位点太少，是无法用于训练一个合适的模型的，这对于全基因组来说是没任何问题的，但外显子区域加起来也就差不多50Mb左右，长度不大，单个样本里面包含的变异数目大约30K-40K。这些位点本来就不多，它们和已知高质量变异集Overlap的就更少了,最终就导致达不到模型训练的最低要求。所以单个样本的WES（或者样本数量较少的WES）都无法使用VQSR进行质控，小Panel的测序数据也是同理。

但随着样本数目的增加，群体中会有更多的变异也在这些外显子区域中被发现，从而增大了这个可用的训练集合，直到满足了最低训练要求，按照经验，通常是30个样本（随着捕获区域的差别，会略有差异），这也是为什么对于WES数据而言，GATK会提到至少需要30个样本才能进行VQSR的原因。

## 3. 样本量重要，还是测序深度重要？

我认为是样本量远比测序深度重要。只要有足够多的样本，我们甚至可以用很低的测序深度（比如1x）获得这些样本中每个人准确的genotype和群体的遗传频谱。这是为什么？

其中一个核心原因是人类这个物种具有单一祖先起源，这也是一个重要的前提假设。但同时我想强调一点，这里的“单一”并不是特指只有一个个体，而是指形成这个群体（比如说现代人，甚至就只是中国的汉族人）的祖先归结起来只有为数不多的若干个部落。在这种情况下，人群多样性的源头实际上就主要来自这些部落之间的基因交流和融合。至于什么是基因交流，大家可以自行脑补。

![]()

另一个核心原因是时间不够。人类其实是一个很年轻的群体，特别是现代智人（我们这一波），遗传的分化历史很短，按照目前估算大约是10万年前才开始。而群体出现遗传差异的动力主要有两个：（1）基因组自身的突变和重组；（2）生殖细胞在形成配子过程中发生的重组。但基因组突变和重组的速率都是很低的，大概只有10^-8次方左右。也就是说一个人因为突变所带来的遗传差异，积累起来大约是30-100个。这个只是序列上的突变（主要是点突变），重组虽然有所不同——它是大范围序列的交换，影响的范围很大，但是一般不认为它直接带来序列突变。我们可以理解为它带来的是突变在整个群体中的扩散和分配。

然而，10万年的时间，差不多只有5000代人，这个数字放在物种遗传的历史上是很短暂的一瞬，这个时间跨度不足以引起整个群体的多样性爆发。对于东亚人来说则更少，目前发表过的研究表明，东亚人的历史更短，大概起源于6万年前，所以你会在千人基因组项目中看到东亚人（特别是汉族人）内部的分化差异极小。最终归结起来，人类这个群体中单倍体的组合数目是非常有限的。

所以如果要揭示一个特定群体的遗传图谱，我们大可不必对全体样本都进行高深度测序，只需要把其中一部分人进行深测获得较高质量的变异集合，然后其他样本则直接使用低深度测序（甚至是定制的芯片测序，不过我更偏向于选择低深度全基因组测序），再结合连锁不平衡遗传定律，我们就完全有能力推断那些没被充分覆盖的区域中的具体基因型，千人基因组和冰岛人就是这样的一个例子。

> GATK的HaplotypeCaller算法实际上也是利用这样的原理实现了更加准确的变异检测的。在变异检测时，GATK会利用所有样本的数据，预先构造出这个群体的Haplotype组合（这应该也是HaplotypeCaller这个名字的由来），以及这个组合中各个单体型在群体中的后验概率，然后再依据每个样本自己的比对数据，通过贝叶斯原理计算出各个样本在每个位点上的基因型和各自基因型的后验概率。如果参与分析的样本足够多，那么理论上它就能够构建出更加准确的Haplotype组合，然后反过来就会提升各个样本的变异检测结果。

## 4. 怎么通过LD衰减距离去看群体的一个遗传多样性呢？

LD本身反应的是一个物种基因组上发生过的重组情况。基因组的重组在每一代都会发生，如果一个群体越古老，那么可以预期它基因组中发生过重组的次数就越多，那么相应的它的LD长度就会越短，从而这个族群的遗传多样性就越高。比如在现代人类中，遗传多样性最高的是非洲人，他们历史最久远，而我们东亚黄种人，多样性则是最低的。如果我们要通过基因芯片对非洲人的某些特征进行全基因组关联分析，那么理论上适合这个群体的芯片密度要比我们黄种人的高。

## 5. 生物信息工程师可以分为多少种类型？

总的来说包含三个大的分类导向：

第一类，**技术导向**，目标是开发更好的算法，思考如何利用数理和计算机等方面的知识提供更好的工具和平台。帮助解决组学问题，比如编写比对算法、组装算法、变异检测算法、质控程序等，当然也包括编写生产级别的数据分析流程（如标准化WGS流程），这一类型的生信工程师解决的是生产工具的问题。

第二类，**数据导向/问题导向**，或者叫“业务”导向——这里的业务包括科学研究和商业应用。主要是解决生物和组学问题、遗传咨询等，如癌症研究、群体遗传学等。这类人更多的是工具的使用者，他们会根据具体的“业务”需要组合最合适的算法和工具来解决问题，这一类人需要较深的生物和基因遗传学知识背景。同时，必须对自己所在的领域有一个完整的认识，知道在什么场景下需要什么数据，应用什么算法，使用什么数理知识和什么工具，才能更好地解决问题——其实这一类人也是真正知道该做什么分析流程的人。

关于这一类生信工程师，或者应该称为“基因组学专家”更加合适，他们包含很多方面，比如群体遗传学、动植物基因组学、进化、肿瘤研究、医学基因检测、消费级基因检测、遗传咨询等。他/她们通常是依据“业务”目标，运用相应的技术手段和工具（包括WGS、WES、RNAseq、甲基化测序、相关组学分析方法等）解决达成目标道路上的问题。这里每一个都可以再进一步展开，总的来说，这个类型是工具的使用方，具体组学问题的解决者。

上面这两类看起来各有特点，掌握的知识点各有侧重，但其实并不能割裂，真正做得好的人，都是两类通吃的（可能只是两强相较，某一类更突出）。因为能深刻理解生物问题和组学问题的人，才能创造出真正合适的工具和流程。

第三类，**资源和人导向**，或者叫“Boss”/PI导向。这些人由于各自成长经历的不同，可能已经和上面的情况有所出入了（很难说会全都懂），他们中有些可能更擅长于去找资源，搭桥，做连接。他们更多的不是解决具体问题，而是尽可能地提出好问题，发现好方向，并为提供解决这些问题创造环境和条件。这一类人其实往往也是第一类和第二类人发展在后面的一个方向。

***

技术交流圈往期精华

* [RNA-Seq是否可以替代WES完成外显子的变异检测?二代测序的四种Read重复是如何产生的? ](https://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798583&idx=1&sn=c8bf9ef0dce441882cd3b80243597756&chksm=83c1d5abb4b65cbdc8b76be8381f8bb6cf5a40b2b83455f8cc7f615c66cd9092c2c5ff2d35ec&scene=21#wechat_redirect)  
* [RNA-seq原始数据质控后，是否要合并PE和SE的比对结果](http://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798518&idx=1&sn=29dfbb2279202ffc353c09916f36c9b8&chksm=83c1d56ab4b65c7c1e4da54edfda880cfc86be8f6d3a8852a8a145b72b45144fabebf0a8a12c&scene=21#wechat_redirect)

* [我是解螺旋的矿工，我热爱生命科学](http://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798476&idx=1&sn=52cd7f2d86d77410d538ab9edbec8df7&chksm=83c1d550b4b65c468425dfc367c23c4da80ebfe7060df2bd314354a3b4690c9425640153ba51&scene=21#wechat_redirect)

* [该如何自学入门生物信息学](http://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798366&idx=1&sn=b545fcea7f82839fa87e9d9e472d1e72&chksm=83c1d4c2b4b65dd4843250c307969ada96c4039f4f528c034620d25b78d8beba2f9cf924bb8a&scene=21#wechat_redirect)

***

本文首发于我的个人公众号：**helixminer（碱基矿工）**

![helixminer-QRCode](https://static.fungenomics.com/images/2021/03/helixminer-mid-red-20210327230714799-20210327230740141.png)

***

这是知识星球：『解螺旋技术交流圈』，是一个我与读者朋友们的私人朋友圈。我有9年前沿而完整的生物信息学、NGS领域的工作经历，在该领域发有多篇Nature级别的科学文章，我也希望借助这个知识星球把自己的一些微薄经验分享给更多对组学感兴趣的伙伴们。

自从星球正式运行以来，已经过去了6个月，星球的成员也已经超过220人了。所分享的主题超过了500个，回答的问题超过了140个，精华70个。我在知识星球上留下的文字估计也已经超过10万字，加上大家的就更多了，相信接下来星球的内容一定还会不断丰富。另外，上周获得了知识星球官方评选的“最优质星球”优秀奖。

这是知识星球上 **第一个真正与基因组学和生物信息学强相关的圈子。**我希望能够借此营造一个高质量的组学知识圈和人脉圈，通过提问、彼此分享、交流经验、心得等，彼此更好地学习生信知识，提升基因组数据分析和解读的能力。

在这里你可以结识到全国优秀的基因组学和生物信息学专家，同时可以分享你的经验、见解和思考，有问题也可以向我提问和圈里的星友们提问。

知识星球邀请链接：[「解螺旋技术交流圈」](https://wx.zsxq.com/mweb/views/joingroup/join_group.html?group_id=518881585444&secret=vcdvs4rdpst7stq4wcvqmlwvogc0ssbn&user_id=28821152428221)

