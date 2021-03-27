---
title: 'GATK4全基因组数据分析最佳实践,我以这篇文章为标志,终结当前WGS系列数据分析的流程主体问题 | 完全代码'
date: 2018-05-11 01:00:00+0800
image: http://image.fungenomics.com/human_genome_project.gif
categories:
    - 生物信息
    - 基因组学
tags:
    - NGS
    - WGS
    - 流程
---

![](http://image.fungenomics.com/human_genome_project.gif)

这是我根据之前的WGS系列和GATK4实践文章进行重新梳理之后确定下来的分析流程，这是一个WGS的最佳实践，它基于GATK4和我的实际经验，稍作修改即可应用到实际的项目中。我以这篇文章为标记，终结当前WGS系列数据分析的主体流程问题。

首先，要提醒大家的是，我在这里所呈现的仅仅只是一个shell流程，它并不符合正规的IT设计规范，你也需要对其参数做适当的修改，比如修改相关软件路径和数据路径，放在这里主要有以下两个目的：

* 第一、我想清楚地告知大家一个完整的WGS数据分析流程里面到底都有啥，执行顺序是怎么样的，先给出一个直接可用的版本供参考，大家在有需要的时候，就不用从零开始了，有据可依，少走一些弯路，可以结合实际项目和流程中的注释信息稍作修改即可使用，提高效率；
* 第二、打下一个最基本的WGS流程基准，为以后在这个基础之上跌代修改提供依据。

另外，GATK4虽然已经正式发布了一段时间，用法也确实和GATK3有些不同，但是差距很小，主要是参数形式和调用形式上的改变（我之前的两篇GATK4实践文章对此已有过演示）。它主要是整合了多个现有工具、新增一部分新功能（包括Spark功能），但基本的核心内容并没有改变。

之前我听说有些同学觉得改换了GATK4之后，担心原来基于GATK3的WGS分析方法就不适用了。这其实有些多虑了，WGS数据分析和解读是我们的目的，GATK是帮助我们达成该目的的一个重要工具，但如果有更好/更合适的工具，我们随时准备替换，甚至重写。所以GATK也好，BWA也罢，对于我们而言都只是“术”，重要的是，我们要知道该如何对数据进行分析和解读，这是根本之“道”。这也是我在写作过程中努力坚持的一个宗旨和原则。

> 另外，我在以下流程的注释中留下了很多重要的信息，以及一些步骤的使用条件，例如某些步骤（比如HaplotypeCaller）可以有多种不同的实现路径，这个可以根据实际的需要进行选择。
> 

好了，现在进入正篇。

我之前已经用四五万字写了一系列的WGS文章，对其中的很多原理和原则都做了详细的解释，这里就不再赘述了，如有需要可以在本文下方的推荐阅读中查阅。另外，在这个流程中我做了一些默认设置，比如，参考序列和GATK4所需的bundle数据都是选择了目前最新的hg38。

那么，接下来，就直接上可用的流程代码吧，我设计为两种模式：

## 第一：单样本模式

你可以把下面这一段代码，写到一个shell脚本中（注意代码中的“\”是接下一行的符号，它后面不能有任何空格或者其它字符），我们这里把文件名定为，wgs_single.sh。这个流程假设你只有一个样本，这个样本只有一对Illumina测序的PE fastq数据文件。流程有6个参数，你可以在代码中清楚地看到，分别是：

* read1的路径
* read2的路径
* Read Group ID
* 测序文库编号
* 样本编号
* 输出目录路径

用起来也比较简单，直接在命令行中执行，并接上以上对应的参数即可，比如：

```bash
$ sh wgs_single.sh read.1.fq.gz read.2.fq.gz Test_RG Test_lib Test_sample Test_outdir
```

以下是完整的流程代码：

![](http://image.fungenomics.com/wgs_single_1.png)
![](http://image.fungenomics.com/wgs_single_2.png)
![](http://image.fungenomics.com/wgs_single_3.png)
![](http://image.fungenomics.com/wgs_single_4.png)

## 第二：多样本模式

以上单样本的流程比较简单直接，如果你没有做成一键式产品的需要，那么基本上是可以用上面的流程一步步完成你的WGS分析的。但随着目前测序行业的发展，大规模人群的测序会变得越来越普遍，并且基于群体的数据分析要远优于单样本，意义也更加深远，因此这种多样本的模式将成为常态。

多样本的流程和单样本相比会有些不同，首先是不适合把执行过程都封装在同一个shell脚本中。在这里我提供了一种实现方式，可供参考，我分为两步：

* 第一步，单独为每个样本生成后续分析所需的中间文件——gVCF文件。这一步中包含了对原始fastq数据的质控、比对、排序、标记重复序列、BQSR和HaplotypeCaller gVCF等过程。这些过程全部都适合在单样本维度下独立完成。值得注意的是，与单样本模式不同，该模式中每个样本的gVCF应该成为这类流程的标配，在后续的步骤中我们可以通过gVCF很方便地完成群体的Joint Calling；

* 第二步，依据第一步完成的gVCF对这个群体进行Joint Calling，从而得到这个群体的变异结果和每个人准确的基因型（Genotype），最后使用VQSR完成变异的质控。这两个步骤其实还包含了许多细节，具体可见我在流程中的注释。

以下是第一步的代码，参数和单样本模式一样。这个代码我们可以把它存放在一个名为wgs_fastq_to_gvcf.sh的文件中：

![](http://image.fungenomics.com/wgs_pop_fastq_to_gvcf_1.png)
![](http://image.fungenomics.com/wgs_pop_fastq_to_gvcf_2.png)
![](http://image.fungenomics.com/wgs_pop_fastq_to_gvcf_3.png)

接着，这是第二步的代码，我们把它存放在一个名为wgs_gvcf_to_vcf.sh的文件中。需要强调的是，它的输入和输出参数需要与第一步保持一致，流程中我对此做了额外的注释：

![](http://image.fungenomics.com/wgs_pop_gvcf_to_vcf_1.png)
![](http://image.fungenomics.com/wgs_pop_gvcf_to_vcf_2.png)
![](http://image.fungenomics.com/wgs_pop_gvcf_to_vcf_3.png)

使用方式可以参考上文“单样本模式”的例子，直接在命令行中完成即可，不再赘述。

![](http://image.fungenomics.com/%E5%88%86%E4%BA%AB.jpg)

## 变异注释

变异的注释这一个步骤，我把它单独拎出来，原因是它完全可以独立于以上的流程。任何SNP和Indel数据（VCF格式），只要有需要你都可以随时完成这个注释。由于比较简单，流程的细节我在这里就不再多说了，最难的其实只是如何安装好VEP和它需要的相关数据集（cachedir目录下的数据）。

```bash
## 使用VEP完成变异的注释
VEP=/your_path_to/ensembl-vep/vep
time $VEP --fasta $reference/Homo_sapiens_assembly38.fasta \
 --vcf --merged --fork 10 --hgvs --force_overwrite --everything \
   --offline --dir_cache /your_path_to/ensembl-vep/cachedir \
   -i $outdir/gatk/${sample}.HC.VQSR.vcf.gz \
   -o $outdir/gatk/${sample}.HC.VQSR.VEP.vcf.gz
```

## 小结

好了，这篇文章就到此结束了，祝你数据分析顺利~


***

## 推荐阅读

*   [GATK4.0和全基因组数据分析实践（上）](http://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798425&idx=1&sn=ae355ed362848578e5c853413f23dfd7&chksm=83c1d505b4b65c13124c9acd210356c4364ec9f5498bbd16fa4475be29811213abb64ea9720f&scene=21#wechat_redirect)
*   [GATK4.0和全基因组数据分析实践（下）](http://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798455&idx=1&sn=67a7407980a57ce138948eb46992b603&chksm=83c1d52bb4b65c3dde31df94e9686654bf616166c7311b531213ebf0010f67a32ce827e677b1&scene=21#wechat_redirect)
*   [从零开始完整学习全基因组测序数据分析：第4节 构建WGS主流程](https://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798296&idx=1&sn=790d0141eec792b25083c63e87fee14c&chksm=83c1d484b4b65d921fd0f17b24e22e17ba76b7e1ca338712298af8bd7532025367d9f47cf630&scene=21#wechat_redirect)
*   [该如何自学入门生物信息学](http://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798366&idx=1&sn=b545fcea7f82839fa87e9d9e472d1e72&chksm=83c1d4c2b4b65dd4843250c307969ada96c4039f4f528c034620d25b78d8beba2f9cf924bb8a&scene=21#wechat_redirect)

***

欢迎关注我的个人公众号：**helixminer（碱基矿工）**

![helixminer-QRCode](https://static.fungenomics.com/images/2021/03/helixminer-mid-red.png)

***

这是我的知识星球：**解螺旋技术交流圈**，是一个与读者朋友们的私人朋友圈，欢迎你的加入。我有9年前沿而完整的生物信息学、NGS领域的工作经历，在该领域发有多篇Nature级别的科学文章。

这是知识星球上 **第一个真正与基因组学和生物信息学强相关的圈子**。我旨在营造一个高质量的组学知识圈和人脉圈，通过提问、彼此分享、交流经验、心得等，彼此更好地学习生信知识，提升基因组数据分析和解读的能力。

在这里你可以结识到全国优秀的基因组学和生物信息学专家，同时可以分享你的经验、见解和思考，有问题也可以向我提问和圈里的星友们提问。

知识星球邀请链接：[「解螺旋技术交流圈」](https://wx.zsxq.com/mweb/views/joingroup/join_group.html?group_id=518881585444&secret=vcdvs4rdpst7stq4wcvqmlwvogc0ssbn&user_id=28821152428221)


