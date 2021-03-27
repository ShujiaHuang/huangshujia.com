---
title: '关于Indel，我该往左还是该往右'
date: 2018-03-01 01:00:00+0800
image: https://static.fungenomics.com/images/2021/03/Left_or_right.jpg
categories:
    - 生物信息
    - 基因组学
tags:
    - NGS
    - 变异检测


---



实质上，两种做法都可以。但是在NGS数据分析中，约定俗成的做法是：左移（left alignment）！然而也有一种情况是例外，当变异使用的是HGVS规则命名的时候，它们却是右移（right alignment）！

在探讨这个问题之前，我们来理解一下为什么Indel需要左移或者右移。

![几米绘本](https://static.fungenomics.com/images/2021/03/jimi-20210327225825714.png)

首先，**会碰到这个问题的场景是我们需要比较多个不同VCF数据集中Indel结果是否一致的时候。**比如，需要把自己项目的VCF数据和标准集比较——这个标准集可以是1000 Genomes Variants（该公开数据集已经左移），也可以是GIAB（Genome In a Bottle）的数据集，也可以是医学检测报证做质评的标准集等等。这个时候我们需要统一Indel的坐标（特别是Deletion的坐标），这是需要「移」目的。

然后，我们来理解一下为什么会出现「移」这个问题的原因。主要有两个：

* 第一，宏观上，基因组本身存在重复性序列，比如在人类基因组中重复性序列占比约50%；
* 第二，微观上，基因组的局部区域存在很多相似序列和短串联重复序列。

在这些地方——特别是原因二，序列比对往往会懵圈，不知道应该在哪个位置上开GAP——Indel（主要是Deletion）在比对上所体现出来的现象就是GAP，才合适。由此，导致的结果就是比对时开GAP的位置每！一！次！都！会！有！差！异！也就是说即使是同一个数据，多次比对，结果也不同。

比如下图，尽管这个区域中实际上只有一个Deletion，但是三次比对的结果，Deletion的断点都不一样，但是从比对结果看，它们的罚分却一模一样，根本无法区分！

![比对懵圈](https://static.fungenomics.com/images/2021/03/repeat_mapping-20210327225825743.jpeg)

理解了上面这一点之后，我们就知道不能在不同数据集中随便对Indel进行直接比较了。我们需要先对这些变异进行规范化的处理（一般指左移），确保所有的Indel的断点都是开在同一个方向上的，用NGS的术语，管这个做法叫：Variant Normalization。

那么，有什么工具可以干这个事情吗？有很多。并且进行左移的合适时机有两个，都有不同的工具可以用。

第一个，在call变异之前，对比对文件SAM/BAM/CRAM直接进行左移，能用的工具很多，比如GATK4.0中的LeftAlignIndels模块。

第二个，已经有了VCF变异集，这个时候GATK中就没有合适的模块来处理了，推荐使用[vt](https://github.com/atks/vt)——一个很不错的变异处理工具，它有一个vt normalize模块，直接就可以对VCF中的Indel断点进行左移。

```
description : normalizes variants in a VCF file.

usage : vt normalize [options] <in.vcf>

options : -o output VCF file [-]
     -d debug [false]
     -q do not print options and summary [false]
     -m warns but does not exit when REF is inconsistent
       with masked reference sequence for non SNPs.
       This overides the -n option [false]
     -n warns but does not exit when REF is inconsistent
       with reference sequence for non SNPs [false]
     -f filter expression []
     -w window size for local sorting of variants [10000]
     -I file containing list of intervals []
     -i intervals []
     -r reference sequence fasta file []
     -? displays help
```

那。。。如果要右移呢？这是一个难得一见的操作，github有一个程序可以用，[这里](https://github.com/counsyl/hgvs/blob/master/pyhgvs/variants.py)。

或许你也好奇，既然右移很少见，那么，为什么会有HGVS这个例外呢？

我也不知道真正的原因，有一种说法是，搞NGS的人和搞医学检测的人对变异的理解有差异，导致各自的处理方式不同，且没有相互沟通过，等发现差异的时候已经迟了。但我不知道这个说法是否属实，姑且听之任之吧。

------------

本文首发于我的个人公众号：**helixminer（碱基矿工）**

![helixminer-QRCode](https://static.fungenomics.com/images/2021/03/helixminer-mid-red-20210327225805055-20210327225826042.png)