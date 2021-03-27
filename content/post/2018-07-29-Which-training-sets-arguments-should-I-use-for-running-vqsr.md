---
title: '我应该如何正确设置GATK VQSR的模型训练参数 |《解螺旋技术交流圈》精华第4期'
date: 2018-07-29 01:00:00+0800
image: http://image.fungenomics.com/bird-s-eye-view-boats-colors.jpg
categories:
    - 生物信息
    - 基因组学
    - 《解螺旋技术交流圈》精华
tags:
    - GATK
---

![](http://image.fungenomics.com/bird-s-eye-view-boats-colors.jpg)

变异的质控，是我们在得到变异数据之后，接下来最重要的一个步骤。通常我们都是使用GATK VQSR模块来完成这个事情，关于VQSR的基本原理我在[这篇文章](https://mp.weixin.qq.com/s/HeIhMeA6GNboQQl1b79arg)中有写，但暂时不算详细。下面是大家经常都会用到的VQSR基本命令（以GATK4为例）：

```bash
## 首先是SNP mode
time $gatk VariantRecalibrator \
   -R $reference/Homo_sapiens_assembly38.fasta \
   -V $outdir/poplation/${outname}.HC.vcf.gz \
   -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $GATK_bundle/hapmap_3.3.hg38.vcf \
   -resource:omini,known=false,training=true,truth=false,prior=12.0 $GATK_bundle/1000G_omni2.5.hg38.vcf \
   -resource:1000G,known=false,training=true,truth=false,prior=10.0 $GATK_bundle/1000G_phase1.snps.high_confidence.hg38.vcf \
   -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $GATK_bundle/dbsnp_146.hg38.vcf \
   -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
   -mode SNP \
   -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 95.0 -tranche 90.0 \
   -rscriptFile $outdir/poplation/${outname}.HC.snps.plots.R \
   --tranches-file $outdir/poplation/${outname}.HC.snps.tranches \
   -O $outdir/poplation/${outname}.HC.snps.recal && \
time $gatk ApplyVQSR \
   -R $reference/Homo_sapiens_assembly38.fasta \
   -V $outdir/poplation/${outname}.HC.vcf.gz \
   --ts_filter_level 99.0 \
   --tranches-file $outdir/poplation/${outname}.HC.snps.tranches \
   -recalFile $outdir/poplation/${outname}.HC.snps.recal \
   -mode SNP \
   -O $outdir/poplation/${outname}.HC.snps.VQSR.vcf.gz && echo "** SNPs VQSR done **"

## 然后是Indel mode
time $gatk VariantRecalibrator \
   -R $reference/Homo_sapiens_assembly38.fasta \
   -input $outdir/poplation/${outname}.HC.snps.VQSR.vcf.gz \
   -resource:mills,known=true,training=true,truth=true,prior=12.0 $GATK_bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf \
   -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $GATK_bundle/dbsnp_146.hg38.vcf \
   -an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
   -mode INDEL \
   --max-gaussians 6 \
   -rscriptFile $outdir/poplation/${outname}.HC.snps.indels.plots.R \
   --tranches-file $outdir/poplation/${outname}.HC.snps.indels.tranches \
   -O $outdir/poplation/${outname}.HC.snps.indels.recal && \
time $gatk ApplyVQSR \
   -R $reference/Homo_sapiens_assembly38.fasta \
   -input $outdir/poplation/${outname}.HC.snps.VQSR.vcf.gz \
   --ts_filter_level 99.0 \
   --tranches-file $outdir/poplation/${outname}.HC.snps.indels.tranches \
   -recalFile $outdir/poplation/${outname}.HC.snps.indels.recal \
   -mode INDEL \
   -O $outdir/poplation/${outname}.HC.VQSR.vcf.gz && echo "** SNPs and Indels VQSR (${outname}.HC.VQSR.vcf.gz finish) done **"
```

在WGS系列文章中大家也可以到类似的程序操作命令。**但是大多数初学者可能并不完全理解GATK VQSR中训练集参数（-resource）的内在含义**，我以前在文章里也缺少对此进行过深入讨论。前段时间在知识星球里有几位小伙伴反复问到了这个问题，我一一作了回答，最后进行了总结。虽然这种参数是GATK VQSR模块所特有，但如果能理解好这些参数以及它们背后包含的意义，**应该有助于我们更好地理解基因组的变异质控，更恰当地使用GATK。**另外，对于生信算法开发者来说，还可以从这样的策略中或多或少地得到一些启发。

所以，这篇文章里，就把我在知识星球中有关这一块的内容整理出来分享给大家。

我们知道，**SNP和Indel的特征是不同的，评价它们的好坏需要使用不同的标准。**因此在使用VQSR进行变异质控的时候，**它们各自的评估模型就需要分开训练和计算**——我想大家在使用GATK分析WGS数据的时候就应该知道了（GATK通过参数-mode SNP或者-mode INDEL来有目的地选择SNP或者Indel，不需要自己去把VCF的SNP和Indel分出来），它们用于训练的数据集也是不完全一样的，所以下面我也按照这两个方面分开总结。

## SNP

对于SNP来说，VQSR的训练集数据目前主要有四个：HapMap，Omni，1000G和dbSNP，参数一般如下：

```bash
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf \
-resource:omni,known=false,training=true,truth=false,prior=12.0 1000G_omni2.5.hg38.vcf \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 1000G_phase1.snps.high_confidence.hg38.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_146.hg38.vcf \
-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
-mode SNP \
```

这个参数的格式，从左到右，按照逗号（,）隔开，分别是：

训练集名字，这个名字是可以随便改动的，但是为了便于交流，一般还是默认按照数据集的名字来设置（如上面的例子）；

* known：该数据是否作为已知变异数据集，用于对变异数据的标注；
* training：该数据是否作为模型训练的数据集，用于训练VQSR模型；
* truth：该数据是否作为验证模型训练的真集数据，这个数据同时还是VQSR训练bad model时自动进行参数选择的重要数据；
* prior：该数据集在VQSR模型训练中的权重，或者叫Prior likelihood（这里转化为Phred-scale，比如20代表的是0.99）。

一共有4个。**但实际上，VQSR训练集的数据可以用-resource参数继续往下添加，有多少就加多少，对于人类WGS/WES数据来说，目前用的主要还是上面的这4个，**但它们其实还是有些旧了（用了8年左右了），现在大人群的数据越来越多，这些数据集估计不久之后也会有变化，比如有些研究也开始加入gnomAD的数据。另外，如果你发现对你来说，这些数据集都不合适，比如都是非人的物种，那么应该自己按照需要构造这些数据集，不一定同样需要四个，只要保证有专门用来训练（training）和验证模型的真集（truth）数据就可以。

在上面的例子里，大家可能也注意到了，对于不同的数据集，与之对应的参数也是不同的，比如hapmap和1000G就不一样，除了prior值不同之外，known、training和truth的参数也不一样。**但这些差异的实质和所代表的意义什么呢？在使用的时候到底是基于什么原则设定了这些参数？**

**毫无疑问的是这些差异的实质一定源自于不同的数据集本身**，所以在逐一解释之前，我想先跟大家一起搞明白HapMap，Omni，1000G和dbSNP这四个数据集分别是什么，它们是怎么来的，以及它们各自都有什么特点。**只要能够弄清楚这几点，那么参数的设置也会不言自明。**

**第一个是HapMap**，它来自国际人类单倍体型图计划，HapMap的名字也是源自于此。这个项目刚启动之时，只有270样本，其中有60个家系。项目一共有三期，到第三期HapMap3的时候这个数据已经扩增到1301个样本了，其中有部分样本和千人基因组项目有重叠。**由于这个数据集包含了大量家系数据，并且有非常严格的质控和严密的实验验证，因此它的准确性是目前公认最高的。**

![](http://image.fungenomics.com/hapmap.jpg)

所以VQSR进行质控模型训练的时候，会将其作为一个很重要的训练集（training=true）。它的权重也会被设置得很高，比如在WGS数据分析中常常设为prior=15——这里的Prior是Prior likelihood的Phred-scale，我们如果把15转换为likelihood，那么就是0.96838。此外，由于它的高准确性，通常还将作为模型验证的一个真集数据（truth=true）。

**第二个是Omni，**这个数据源自Illumina的Omni基因型芯片，大概2.5百万个位点（我在知识星球中说是1M，这里纠正一下，应该是2.5M），它的验证结果常常作为基因型的金标准。比如用Omni芯片对千人基因组数据进行了验证的结果：1000G_omni2.5.hg38.vcf（这里hg38是参考序列版本），它也是一个高可信的变异结果，我们在VQSR模型训练的时候，同样可以为其设置很高的权重，一般为Prior=12(likelihood为0.9369)。**通常情况下也可以把它作为验证结果的真集数据，**但我这里所举的例子有些保守，把它设置为truth=false了，大家不必效仿，**假如你没太大的“洁癖”把它设置为true都是没问题的（这也是GATK最佳实践的一般做法）。**

**第三个是1000G，**这个数据从名字看我们也很熟悉，它就是千人基因组计划（1000 genomes project）质控后的变异数据，目前也是第三期，一共包含了2504个人的数据。通常来说质控后，它包含的绝大部分都是真实的变异，但由于没办法做全面的实验验证，并不能排除含有少部分假阳的结果。所以模型训练时给的权重虽然比较高——prior=10(likelihood为0.9)，但是一般就不作为模型验证的真集数据了，即truth=false。

﻿﻿![](http://image.fungenomics.com/1kgp_nature.jpeg)

**第四个是dbSNP。**说到dbSNP，**这是一个绝对不可以作为训练集位点的数据——太脏了**，为什么这么说呢？因为，**dbSNP收集的数据，实际都是研究者们发表了相关文章提交上来的变异，这些变异很多是没做过严格验证的，很多甚至还是假的，在没被反复验证之前，是不可信的。**因此，不会把它们作为模型的训练集，更不会把它作为真集看待（training=false,truth=false），权重也一般设置得很低，比如这里是prior=2（差不多才0.37）。dbSNP的唯一作用就是用于标注我们的变异集中哪些是已经在其它研究中出现过的——即属于已经被发现过的（已知）变异，给这些已知的变异位点标上RS id。

## Indel

对于Indel来说，VQSR模型的训练参数只有两个：Mill和dbSNP。

```bash
-resource:mills,known=true,training=true,truth=true,prior=12.0 Mills_and_1000G_gold_standard.indels.hg38.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 dbsnp_146.hg38.vcf \
-an DP -an QD -an FS -an SOR -an ReadPosRankSum -an MQRankSum \
-mode INDEL \
```

**第一个是Mills，**对于Indel来说能正在算得上真集的并不多，Mills_and_1000G_gold_standard.indels.hg38.vcf算是其中一个，并被专门做过验证。但其实Indel并不那么容易验证，很多时候也是挑一些比较容易验证的结果。但不管如何，这是目前最佳的一个！所以权重各方面也都设置的比较高，比如prior=12，并且truth=true。

**第二个还是dbSNP，**这个就和上面SNP模式下的作用是一样的。不过假如这一步对Indel进行VQSR的VCF数据是顺着上面SNP VQSR后下来的话，那么这个dbSNP的参数可以省略，因为已知变异的标注已经在SNP model下做好了。我在[WGS系列的这篇文章](http://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798570&idx=1&sn=9445ed66605c8c62595458fa6b561e1d&chksm=83c1d5b6b4b65ca0e3ec97082f65dbdcb6948aaae7560353115a5aa2dbd6ed7f8e233f34c3f1&scene=21#wechat_redirect)里，给出的流程就是这种情况。

## 一些需要额外注意的地方

GATK VQSR在执行的时候要基于全基因组的所有变异数据，而不能拆分不同的染色体分别取执行，否则会导致各个染色体所训练的质控模型不一致。另外，除了要区分SNP和Indel模式之外，GATK VQSR分为两个步骤进行：VariantRecalibrator 和 ApplyVQSR，两者缺一不可。VariantRecalibrator用来进行模型计算，获得数据的情况，ApplyVQSR则是根据我们设定的ts_filter_level参数，最终过滤得到数据，这个参数基于我们对模型真集数据的灵敏度和特异性来确定，一般会设置为99.0%（比如上面例子）。但事实上，这个比例，并不是说它是最合适的，设定为99或者99.5本身并没有什么理论证明，更多还是一种约定俗成，或者大部分情况下，看到设置为这个值，结果看起来是合理的，而且是在正常人样本数据得到的认识。所以，有时候还是需要具体问题具体分析，多看VQSR得到的tranche图，特别是对于非人物种，如果你觉得95.0%也很合适，那么也不一定非得是99%。

***

**你还可以看**

* [人类基因组的Phasing原理是什么？](http://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798506&idx=1&sn=773a7db5dfd002c53fead86625266af5&chksm=83c1d576b4b65c6057f93a00d62ea41364ccdf53429b1820442f5855a104dde68f4325376722&scene=21#wechat_redirect)

* [样本量重要，还是测序深度重要? 生物信息工程师可以分为多少种类型? ](http://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798604&idx=1&sn=48a2e16c6899880fc161cd071fd76651&chksm=83c1d5d0b4b65cc68d27f58c16240a547e7180648a3965afca0917c0f59e4ad42cf62fc4ec46&scene=21#wechat_redirect)

* [GATK4全基因组数据分析最佳实践 ，我以这篇文章为标志，终结当前WGS系列数据分析的流程主体问题](http://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798570&idx=1&sn=9445ed66605c8c62595458fa6b561e1d&chksm=83c1d5b6b4b65ca0e3ec97082f65dbdcb6948aaae7560353115a5aa2dbd6ed7f8e233f34c3f1&scene=21#wechat_redirect)

* [该如何自学入门生物信息学](http://mp.weixin.qq.com/s?__biz=MzAxOTUxOTM0Nw==&mid=2649798366&idx=1&sn=b545fcea7f82839fa87e9d9e472d1e72&chksm=83c1d4c2b4b65dd4843250c307969ada96c4039f4f528c034620d25b78d8beba2f9cf924bb8a&scene=21#wechat_redirect)

***

本文首发于我的个人公众号：**helixminer（碱基矿工）**

![helixminer-QRCode](https://static.fungenomics.com/images/2021/03/helixminer-mid-red.png)

***

这是知识星球：『**解螺旋技术交流圈**』，是一个我与读者朋友们的私人朋友圈。我有9年前沿而完整的生物信息学、NGS领域的工作经历，在该领域发有多篇Nature级别的科学文章，我也希望借助这个知识星球把自己的一些微薄经验分享给更多对组学感兴趣的伙伴们。

这是知识星球上 **第一个真正与基因组学和生物信息学强相关的圈子**。我希望能够借此营造一个高质量的组学知识圈和人脉圈，通过提问、彼此分享、交流经验、心得等，**彼此更好地学习生信知识，提升基因组数据分析和解读的能力**。

在这里你可以结识到全国优秀的基因组学和生物信息学专家，同时可以分享你的经验、见解和思考，有问题也可以向我提问和圈里的星友们提问。

知识星球邀请链接：[「解螺旋技术交流圈」](https://wx.zsxq.com/mweb/views/joingroup/join_group.html?group_id=518881585444&secret=vcdvs4rdpst7stq4wcvqmlwvogc0ssbn&user_id=28821152428221)

