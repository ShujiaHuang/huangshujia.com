---
title: "基于 Cell ID 的单细胞基因指纹特征提取和细胞身份识别的新方法"
description: "跨领域的应用，微创新的作用"
date: "2021-06-26 01:00:00+0800"
image: https://static.fungenomics.com/images/2021/06/pexels-photo-3225618-20210621225816188.jpeg
categories:
    - 单细胞组学
    - 文献阅读
---

我想分享一篇今年四月份发表在 《Nature biotechnology》 上的文章，题目是 “Gene signature extraction and cell identity recognition at the single-cell level with Cell-ID”，翻译出来是 “基于Cell ID的单细胞水平基因指纹特征提取和细胞身份识别方法”。

![image-20210615135757667](https://static.fungenomics.com/images/2021/06/image-20210615135757667.png)

我们知道单细胞RNA测序技术（scRNA-seq）的应用已经越来越广泛。在研究人类器官组织和细胞的类型上，scRNA-seq 是一个很好的技术解决方案。目前比较有代表性的研究项目包括：人类细胞图谱项目、美国国立卫生研究院（NIH）主导的人类生物分子图谱计划和 LifeTime 项目等。

单细胞研究的一个重要目的是揭示细胞之间复杂且丰富的异质性特征。但这个领域一直以来都有一个挑战，那就是 **scRNA-seq 数据的维度和噪声都比较高，这导致对细胞异质性的研究变得十分复杂，这个问题在很大程度上也制约了细胞多样性的研究**。

一般来说降低数据维度可以提高信噪比，也就是用少量但是显著的特征描述细胞。目前在这方面用得最广的是PCA、ICA、tSNE和UMAP。但这些方法本质上都是基于聚类来实现的，它们在计算细胞特征的过程中寻找的是细胞之间的共性差异，**因此这些方法所揭示的异质性实际上停留在细胞亚群这个较粗的颗粒度上，而不是精确到单个细胞的水平**。

但对于细胞转录异质性的详尽研究实际上需要我们对每一个细胞的基因表达特征都进行可靠的统计和分析。所以，**如何能够可靠地表征每一个细胞的基因特征就成为了单细胞转录组学研究中一个尚未被完整解决的科学问题**。

这也是 NBT 这篇文章所要解决的问题。这篇文章的作者来自法国巴黎大学，**他们提出了一个称为 Cell-ID 的无聚类多元统计方法**。这个方法可以从单细胞测序数据中将每个细胞的基因特征有效地提取出来，并且还能够横跨不同的数据集对不同的细胞类型进行注释和匹配，发现未知或罕见的细胞类型和细胞状态。

从原理来说，Cell-ID 背后所依据的方法是统计学中的**多重对应分析法（Multiple correspondence analysis, 简称 MCA）**，这是一种变量统计分析方法，它可以用来分析多个高纬度变量（比如基因表达量）之间的关联以及和多个低纬度观察值（比如细胞）之间的对应关系。

>  MCA方法经常被应用到社会科学领域的研究之中，用来调查测试对象对不同问题的态度一致性，法国和日本对这个方法的使用尤为普遍，可能也是这个原因法国巴黎大学的研究人员率先将这个方法移植到了单细胞基因指纹特征的提取上。

MCA 本质上也是一个针对多个分类变量的降维方法，就如同针对定量变量的主成分分析法一样，最终目的是让同个类别的对象将紧靠在一起，而不同类别的对象远远分开（但不同之处在于MCA除了降维之外，还可以检测多变量之间的关联关系）。

但要注意 MCA 仅适用于从定性变量（也就是“分类变量”）中得出统计结论，**所以在应用 MCA 之前需要先将定量变量转化为分类变量**，例如将连续型变量标准化之后取它们的统计分位数来作为分类变量。当数据集完全表示为分类变量之后，就可以构建相关的数据矩阵进行 MCA 分析了。

文章这里是通过线性变换的方式将各个细胞的基因表达量转化到0和1的范围之内，这样就可以在 MCA 的数学框架之下对细胞的特征进行定量分析了，这一点也是这篇文章中一个较为巧妙的处理方法，下面的图1是 Cell-ID 的原理概述。

![图1. Cell-ID方法概述](https://static.fungenomics.com/images/2021/06/image-20210617084300205.png)

从图1.a可以看出，Cell-ID 通过 MCA 实现了基因表达矩阵的降维，图中细胞（黄色圈圈）和基因（黑色“+”）都投影到一个共同的正交空间中，这个图也叫做 MCA 双标图（MCA biplot）。**在这样的正交空间中，基因离细胞越近，那么就代表它对某个细胞的特异性越高**。因此，可以在 MCA 空间中，将细胞上的基因与该细胞的距离作排序，**排序靠前的基因就可以作为这个细胞的基因特征**，或者称作该细胞的基因指纹，看作是这个细胞的一个独特身份证。同时，每个细胞的基因特征本身也是一个很有价值的数据，所以还可以单独将它们提取出来构成一个单细胞基因特征数据集用于进行下游分析，比如图1b中的功能研究等。

那么关于 Cell-ID 的原理概述就到此为止了，当然具体的数学细节我在这里无法展开，因为该部分的细节对我来说也还有不清楚的地方，还需要做更多的数学演绎才行，当然这其中最重要的就是 MCA 的原理（在很多多元统计学的书本中有该原理的数学描述）。

## 评估 Cell-ID 的有效性

接下来要对 Cell-ID 方法的有效性进行综合评估，**这个评估的方法和结果贯穿全文也是文章的一个重点**。

研究人员首先模拟生成了 100 个 scRNA-seq 数据，然后在这个数据集上分析基于MCA 降维的细胞和基因表达的一致性。

这个一致性评估从三个层面来进行：

1. 通过计算 Spearman 相关性系数的方法，分别评估 MCA 方法和常用的 PCA 方法降维之后前10个主成分的结果相关性（如 Supplementary Fig 1），可以看到各个PC之间的相关性都很高，接近于 1；

![image-20210617202644406](https://static.fungenomics.com/images/2021/06/image-20210617202644406.png)

2. 通过近邻法进行对比。具体来说是对由 MCA 方法所获得的每个细胞中的基因排列和 MCA 空间中相邻的另外50个细胞进行对比。主要是比较他们的基因表达量是否一致，对比的结果也在 Supplymentary Fig1a.b 中展示出来了，可以看出一致性情况也是很好的，**也就是说彼此相邻的细胞，它们的基因特征也相似**；

3. 进一步验证发现基于 MCA 得到的细胞基因特征信息即使是在发生高dropout现象的 scRNA-seq 数据集里也依然有很好的鲁棒性。

> 这里我补充解释一下什么是“Dropout现象”：Dropout 现象是 scRNA-seq 中常发生的一个事件，意思是基因表达信息漏测。原因是很多在表达的基因，由于每个细胞中 mRNA 序列起始量较低或者测序技术的原因而没有被检测到，这部分基因的表达信息就被漏掉了。dropout 现象所导致的数据丢失，会影响下游的数据分析，如何解决这个问题也是单细胞组学所面临的一个挑战。
> https://www.linkresearcher.com/theses/b57bbc38-da8c-463b-8c91-4d56c3101ac4 

除了使用模拟数据之外，研究人员接下来使用两组独立的人血单核细胞对 Cell-ID 的有效性做更进一步的评估，这两组细胞分别是：（1）基于 CITE-seq 方案得到的脐血单核细胞（CBMCs）和（2）通过 REAP-seq 方案得到的外周血单核细胞，这两个方案都是通过检测单个细胞的蛋白标记物水平对单个细胞的特征进行了注释。

>这相当于是通过实验检测的方法，得到了一个单细胞类型特征的参考数据。

通过对比分析之后，可以发现 Cell-ID 的基因指纹特征在对应类型细胞的基因上都有显著的富集情况（图 2.a），这个富集可以说明 Cell-ID 得到的细胞基因特征和真实结果是具有高度一致性。

![图2. Cell-ID通过预先建立的标记列表识别人类CBMCs细胞类型](https://static.fungenomics.com/images/2021/06/image-20210617210451326.png)
图2. Cell-ID通过预先建立的标记列表识别人类CBMCs细胞类型

从具体数字上来说，这两个数据集中 Cell-ID 的识别精确度（Precision）分别达到了 87% 和 90%，召回率（Recall）达到了 84% 和 73%。这个结果想要告诉我们的是，文章所提出的 **Cell-ID 能够很好地提取每个细胞的基因指纹并用来识别不同的细胞类型**。

除此之外，Cell-ID 甚至还能识别正在分化的细胞亚型，例如文章图 2c,d 所示的那样，Cell-ID 捕获到了造血干细胞的分化亚型，而且即使是罕见的细胞状态也可以被 Cell-ID 识别出来。

### 同类细胞的可重复性识别评估

紧接着，研究人员进一步评估了 Cell-ID 对同一组织不同批次 scRNA-seq 数据集中识别相同细胞类型的能力。如文章图3所示（如下），主要分析了来自多个不同供体、不同测序平台所产生的人类胰岛和人类以及小鼠气道上皮细胞的数据集。结果发现，Cell-ID 的整体性能和有效性与已经发表的方法相当。精确度和召回率都很高，其中精确度大于 92%，召回率也高于 75%。

![图3. Cell-ID对同一或不同来源组织、种内和种间的scRNA-seq数据集的细胞匹配表现](https://static.fungenomics.com/images/2021/06/image-20210617221723547.png)
图3. Cell-ID对同一或不同来源组织、种内和种间的scRNA-seq数据集的细胞匹配表现

### 跨组织的细胞类型识别能力评估

然后，评估 Cell-ID 在不同组织来源的 scRNA-seq 数据中识别同一种细胞类型的能力。

文章还是用气道上皮细胞作为例子（图3 展示了这个过程）。基于在气道上皮细胞中获得的无偏基因指纹特征，Cell-ID 识别出了肠上皮中的刷状/族状细胞、内分泌细胞和杯状细胞，而且精度高达90%、召回率达到73%。对比之后发现，这个精度已经优于已经发表过的方法（图3c、d可以查看更加具体的对比结果）。

另外，他们还使用 Cell-ID 对两个独立的嗅上皮细胞数据集做细胞类型的扫描和识别，同时对比了来自气道和肠道上皮的族状细胞特征，**结果还识别出了推测中罕见的、未分类的SCCs细胞，即孤立化学感觉细胞**(如图3e、f所示)。

### 跨测序平台的评估

这是对 Cell-ID 的最后一项评估，评测了它在不同的单细胞组学平台上的表现，同时验证它在基因指纹特征识别上的可重复性（文章图4展示了这个具体的过程）。这个评估所用到的数据主要来自于雄性小鼠细胞图谱中的 scRNA-seq 数据和小鼠 ATAC 图谱中单细胞 ATAC-seq 数据。分析结果也显示，**Cell-ID 对来源于 scRNA-seq 和 ATAC-seq 的数据得出的细胞类型匹配度都很好，F1 分值也都比较高，并且要优于当前已发表的其他方法**。

![图4. Cell-ID对来源于不同单细胞组学技术的独立数据集上细胞间匹配能力的评估](https://static.fungenomics.com/images/2021/06/image-20210617222801065.png)
图4. Cell-ID对来源于不同单细胞组学技术的独立数据集上细胞间匹配能力的评估

所以，综合来说，这篇文章所提出的 **Cell-ID 可以非常量化地提取并注释细胞的基因特征用于表征不同的细胞类型，并且能够在不同的供体、器官组织、物种和单细胞测序平台中得到有效的重复和验证**。这样的一个方法可以改善我们在单细胞水平的生物学方面的研究和解释力，可以更好地发现以前未被表征出来的罕见细胞类型或者细胞状态，而且这个方法还为跨组织、跨生物体的细胞类型研究以及系统多组学研究奠定基础，意义非凡。

## 启发

读了这一篇文章之后，我想最后再谈一点它带给我的一个启示。

实际上，这篇文章所用到的统计学方法并不算很新颖，它只是将一个在其他学科中用得比较广的方法复用到单细胞组学中来，并没有创造一种全新的统计学算法，但是却在单细胞组学领域取得了很好的效果，**可见微创新同样是获得重要科研成果的有效手段，甚至还是一个高效的手段，未必都得好高骛远，还是要因时制宜**，当然他们的先发优势也是这里面一个重要的加分项。

文章提出的 Cell-ID 从原理上来说，我觉得还可以应用到肠道菌群基因组的研究之中，可以用类似于 Cell-ID 的思路设计一个 Meta-ID 来对不同的肠道菌群进行特征表达和识别，这样一来应该也能够进一步提升当前的肠道菌群多组学研究。

*参考文献*

> Cortal A, Martignetti L, Six E, Rausell A. Gene signature extraction and cell identity recognition at the single-cell level with Cell-ID [published online ahead of print, 2021 Apr 29].Nat Biotechnol.2021;10.1038/s41587-021-00896-6.

## 订阅

关注我的个人公众号：**helixminer（碱基矿工）**

![helixminer（碱基矿工）](https://static.fungenomics.com/images/2021/03/helixminer-mid-red.png)
