---
title: '如何使用Shapeit2对人类基因组数据进行phasing'
date: 2018-06-10 01:00:00+0800
description: Shapeit2 教程
image: http://image.fungenomics.com/Van_Gogh_-_Starry_Night.jpg
categories:
    - phasing
tags:
    - shapeit2
---

![](http://image.fungenomics.com/Van_Gogh_-_Starry_Night.jpg)

在上一篇文章中，分享了有关[基因组Phasing的原理](http://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798506&idx=1&sn=773a7db5dfd002c53fead86625266af5&chksm=83c1d576b4b65c6057f93a00d62ea41364ccdf53429b1820442f5855a104dde68f4325376722&scene=21#wechat_redirect)，一共有三种，分别是：家系关系分型(Related individuals Phasing)、群体LD分型(LD Phasing)和物理分型(Physical Phasing)。目前用的比较多的就是基于群体LD的分型，也是我们接触比较多的Phasing分析（另外两种都对数据有一些额外的要求），同时也有比较好用的软件工具，比如很有代表性的Shapeit、beagle和STITCH。在这篇文章里，主要来介绍如何用Shapeit实现基因组的Phasing，另外两个如果有必要的话以后再来补充。

Shapeit是一个专门用于推断基因组单体型（Phasing）的软件，它和beagle一样是当前用得最多的两个基于群体LD进行单倍型推断的软件，使用场景和算法彼此间大同小异。它目前的最新版是Shapeit3，但是常用的还是Shapeit2，也是在千人基因组项目中主要应用的版本。而Shapeit3主要是针对超大规模人群，一般是量级在几万人规模的基因组会更加合适，都是牛津大学的团队开发的，这个版本3可以说是为他们国家的GenomicsEngland计划定制的，这是一个要测10万英国人基因组的大型项目——也是目前世界上推得最快的国家级基因组计划。

![](http://image.fungenomics.com/shapeit2.jpg)

由于主要是说实操，所以这篇文章的内容比较简单，有关Phasing的原理和意义都在[上一篇文章](http://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798506&idx=1&sn=773a7db5dfd002c53fead86625266af5&chksm=83c1d576b4b65c6057f93a00d62ea41364ccdf53429b1820442f5855a104dde68f4325376722&scene=21#wechat_redirect)里仔细讲过了，所以本文篇幅也不长，旨在用具体的例子演示如何使用好这个工具完成基因组数据的Phasing，并构造出Reference panel的过程。

## 首先，准备文件

整个Phasing过程我们需要4个文件：变异数据集（VCF 格式）、样本信息文件(sample.ped)、genetic map和参考序列（fasta格式）。

关于genetic map文件，需要单独拿出来做些说明。这个文件所记录的是基因组中各个位点的重组率和彼此间物理距离的关系，它是一份比较固定的数据文件。目前通用的版本来自于2005年HapMap——人类单体型计划的成果，可以在NCBI上下载，虽然版本的年代比较久远了，但是目前也没有更好的单倍体成果能够去代替它。它是Shapeit2 完成Phasing最重要的一个文件。目前NCBI上[下载的genetic map对应的参考序列版本是b37(即hg19)](ftp://ftp.ncbi.nlm.nih.gov/hapmap/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz)。如果我们的参考序列版本高于hg19，比如现在最新的GRCh38，那么只能进行liftover，把位点坐标转化到最新的参考序列上。

新版本的参考序列会修正旧版本参考序列中的一些问题，包括补齐某些区域的序列或者调整某些区域中序列的顺序。liftover是通过两个不同版本的参考序列进行长序列比对后，获得两个版本间坐标上的映射关系，然后进行坐标转换。既然是通过这种方式获得的，就难免存在一些错误，因此liftover之后，需要对它进行一些过滤处理。主要是在genetic map中那些位点顺序发生交叉的位点，是liftover的错误所导致的，要去掉。

```bash
position COMBINED_rate(cM/Mb) Genetic_Map(cM) 35326 0.251801 0.000000
35411 0.482009 0.000021
40483 0.598191 0.002466
40852 0.599253 0.002687
41421 0.592293 0.003028
41892 0.591345 0.003307
42920 0.622736 0.003915
43259 0.778280 0.004126
44167 1.380848 0.004833
```
以上，是genetic map文件中所看到的内容例子。格式上比较清晰，第一列是位点距离，也就是位点坐标之间的物理距离；第二列是重组率（通过大规模家系数据测量得到）；第三列就是Genetic map的数值，每一个位点都有一个值，计算公式是前一个位点的Genetic Map + physical_distance × recombination_rate。

应该注意到的是，liftover之后，原来genetic map中两个位点之间的重组率（recombination rate）依然是不变的。这其实也很好理解，参考序列之所以需要升级，是因为旧版本的结果并非是百分百符合真实情形的，随着技术的进步，我们会不断去升级逼近真实的序列情况。但重组率是根据大规模家系数据的重组情况来计算的，它是真实情况反映出来的现象，因此即便参考序列版本改变了，它的值也不需要改变。但对于两个位点之间的物理距离（physical distance）来说就不同了，liftover之后，这个距离是有可能发生改变的。

另外一个重要的数据是：样本信息文件，一般称为PED文件，格式如下：

```bash
1009  1009-01  0  0  1
1009  1009-02  0  0  2
1009  1009-06  1009-01  1009-02  2
1030  1030-01  0  0  1
1030  1030-02  0  0  2
1030  1030-06  1030-01  1030-02  1
```
第一列是家系编号，每一个家系有一个唯一编号，可以我们自己人为设置；
第二列是该样本的编号，如果该样本是家系中的小孩，那么需要分别在第三列和第四列中，给出他/她的父亲样本编号和母亲样本编号，而如果该样本是家系数据中的父母，则只需要用0填充即可，如上例子；
第五列是性别信息，一般用1和2分别代表Male和Female；

变异数据（VCF）和参考序列文件就不必多说了。准备好以上文件之后接下来就是主要的流程步骤了。

## 第一步，将vcf转化为bed/bim/fam

bed/bim/fam是Shapeit2完成Phasing所需的三个谱系文件格式。我们需要首先把VCF转为这三个文件，做这一个转换有多个选择，可以用plink，也可以用GATK3（它也有相关模块能够完成这个转换）。我这里就直接用plink来进行转换了：

```bash
plink=/com/extra/testing/bin/plink time  $plink \
  --vcf chr22.vcf \
  --vcf-half-call missing \
  --keep-allele-order \
  --a1-allele chr22.vcf 4 3 '#' \
  --vcf-idspace-to _ \
  --allow-extra-chr 0 \
  --split-x b38 no-fail \
  --make-bed \
  --noweb \
  --out chr22 && echo  "** done **" && sed  's/^chr//g' chr22.bim > t.bim && mv -f t.bim chr22.bim
```

可以看到，我在最后多加了一小步：将原来输出的.bim文件中第一列的chr22换成了22。之所以要费这个小周折，是因为接下来用plink质控的时候，它有点傻，竟然不认识以chr开头的染色体ID！如果没有这个小操作，会碰到`ERROR: Problem reading BIM file, line 1`的报错，从而导致流程中断。

## 第二步，把基因型丢失率（genotype missing rate）较高的位点和含有孟德尔错误的位点过滤掉

基因型丢失率或者孟德尔错误过高的位点，它们的结果可能只是技术噪音(如复杂区域的比对错误，测序PCR导致覆盖偏向)而已，这样会导致基因型推断出错。把这些有问题的位点放入模型，就会导致Phasing的错误率的增高。

```bash
plink=/com/extra/testing/bin/plink time  $plink \
   --noweb \
   --bfile chr22 \
   --keep-allele-order \
   --me 1 1 \
   --set-me-missing \
   --make-bed \
   --out chr22.nomendel &&  echo  "** nomendel done **" && \
time $plink \
   --noweb \
   --bfile chr22.nomendel \
   --keep-allele-order \
   --geno 0.05 \
   --make-bed \
   --out chr22.nomendel.filter && echo "** fileter done **"
```

## 第三步，Phasing

最后阶段了，这里分成两小步来完成：phasing和输出格式转换：

```bash
# phasing  time shapeit2 \
 --duohmm \
 -W 5 \
 --input-bed chr22.nomendel.filter.bed chr22.nomendel.filter.bim chr22.nomendel.filter.fam \
 --input-map genetic_map.chr22.txt \
 -O hapData \
 --thread 1 &&  echo  "** panel  done **"

# 格式转换
time shapeit2 -convert \
 --input-haps hapData \
 --output-vcf chr22.haps.vcf \
 --output-ref chr22.phased.hap chr22.phased.leg chr22.phased.sam && echo "** all done **"
```
大家可能也注意到了，在以上输出的结果中，我设置了两个：--output-vcf和--output-ref。chr22.haps.vcf是Phasing之后的结果，它是一个群体级别的全基因组单体型集合，可以作为我们常说的Reference panel。不过用于Imputation的Reference panel除了这个VCF格式之外，还可以有其它格式，比如这里--output-ref的三个输出文件，也是常用的Reference panel文件。

## 小结

其实从讲Phasing原理的那一刻开始，我的GWAS系列文章就已经在悄悄开始了，接下来当然还有Imputation相关的内容，这两块都只是整个GWAS系列的铺垫。做过GWAS的同学都应该知道，Phasing和Imputation是一定少不了的。

## 致谢

感谢解螺旋的羊在本文写作中的技术讨论和问题解释。

* * *

## 推荐阅读

*   [人类基因组的Phasing原理是什么？](http://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798506&idx=1&sn=773a7db5dfd002c53fead86625266af5&chksm=83c1d576b4b65c6057f93a00d62ea41364ccdf53429b1820442f5855a104dde68f4325376722&scene=21#wechat_redirect)

*   [GATK4.0和全基因组数据分析实践（上）](http://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798425&idx=1&sn=ae355ed362848578e5c853413f23dfd7&chksm=83c1d505b4b65c13124c9acd210356c4364ec9f5498bbd16fa4475be29811213abb64ea9720f&scene=21#wechat_redirect)

*   [GATK4.0和全基因组数据分析实践（下）](http://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798455&idx=1&sn=67a7407980a57ce138948eb46992b603&chksm=83c1d52bb4b65c3dde31df94e9686654bf616166c7311b531213ebf0010f67a32ce827e677b1&scene=21#wechat_redirect)

*   [该如何自学入门生物信息学](http://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798366&idx=1&sn=b545fcea7f82839fa87e9d9e472d1e72&chksm=83c1d4c2b4b65dd4843250c307969ada96c4039f4f528c034620d25b78d8beba2f9cf924bb8a&scene=21#wechat_redirect)

***

欢迎关注我的个人公众号：**helixminer（碱基矿工）**

![helixminer-QRCode](https://static.fungenomics.com/images/2021/03/helixminer-mid-red.png)

***
这是知识星球：『解螺旋技术交流圈』，是一个我与读者朋友们的私人朋友圈。我有9年前沿而完整的生物信息学、NGS领域的工作经历，在该领域发有多篇Nature级别的科学文章，我也希望借助这个知识星球把自己的一些微薄经验分享给更多对组学感兴趣的伙伴们。

自从星球正式运行以来，已经过去了6个月，星球的成员也已经超过220人了。所分享的主题超过了500个，回答的问题超过了140个，精华70个。我在知识星球上留下的文字估计也已经超过10万字，加上大家的就更多了，相信接下来星球的内容一定还会不断丰富。另外，上周获得了知识星球官方评选的“最优质星球”优秀奖。

这是知识星球上 **第一个真正与基因组学和生物信息学强相关的圈子**。我希望能够借此营造一个高质量的组学知识圈和人脉圈，通过提问、彼此分享、交流经验、心得等，彼此更好地学习生信知识，提升基因组数据分析和解读的能力。

在这里你可以结识到全国优秀的基因组学和生物信息学专家，同时可以分享你的经验、见解和思考，有问题也可以向我提问和圈里的星友们提问。

知识星球邀请链接：[「解螺旋技术交流圈」](https://link.jianshu.com/?t=https%3A%2F%2Fwx.zsxq.com%2Fmweb%2Fviews%2Fjoingroup%2Fjoin_group.html%3Fgroup_id%3D518881585444%26secret%3Dvcdvs4rdpst7stq4wcvqmlwvogc0ssbn%26user_id%3D28821152428221)

