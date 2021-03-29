---
title: GATK中如何计算Inbreeding coefficient（近交系数）
date: 2017-05-28 01:00:00+0800
image: https://static.fungenomics.com/images/2021/03/ic_cover.jpg
categories:
    - 生物信息
tags:
    - GATK
---

关于近交系数是什么的定义，除了英文资料，中文上也给出了清晰的定义，这里引用一下：

> 近交系数（inbreeding coefficient）是指根据近亲交配的世代数，将基因的纯化程度用百分数来表示即为近交系数，也指个体由于近交而造成异质基因减少时，同质基因或纯合子所占的百分比也叫近交系数，普遍以F或f来表示。

GATK近交系数的计算程序在github上可以找到：[AS_InbreedingCoeff.java](https://github.com/broadgsa/gatk-protected/blob/f185a75e1c49fb4c039511e61254da0509833ee9/protected/gatk-tools-protected/src/main/java/org/broadinstitute/gatk/tools/walkers/annotator/AS_InbreedingCoeff.java)

代码不短，但计算很简单，我主要说展示一下这个计算的核心部分并在代码中做些注释，如下：

```java
    protected double calculateIC(final VariantContext vc, final Allele altAllele) {
        final int AN = vc.getCalledChrCount();
        final double altAF;

        final double hetCount = heterozygosityUtils.getHetCount(vc, altAllele);

        final double F;
        //shortcut to get a value closer to the non-alleleSpecific value for bialleleics
        if (vc.isBiallelic()) {
            double refAC = heterozygosityUtils.getAlleleCount(vc, vc.getReference());
            double altAC = heterozygosityUtils.getAlleleCount(vc, altAllele);
            double refAF = refAC/(altAC+refAC);
            altAF = 1 - refAF;
            F = 1.0 - (hetCount / (2.0 * refAF * altAF * (double) heterozygosityUtils.getSampleCount())); // inbreeding coefficient
        } else {
            //compare number of hets for this allele (and any other second allele) with the expectation based on AFs
            //derive the altAF from the likelihoods to account for any accumulation of fractional counts from non-primary likelihoods,
            //e.g. for a GQ10 variant, the probability of the call will be ~0.9 and the second best call will be ~0.1 so adding up those 0.1s for het counts can dramatically change the AF compared with integer counts
            altAF = heterozygosityUtils.getAlleleCount(vc, altAllele)/ (double) AN;

            // 计算inbreeding coefficient
            F = 1.0 - (hetCount / (2.0 * (1 - altAF) * altAF * (double) heterozygosityUtils.getSampleCount())); // heterozygosityUtils.getSampleCount() 获取总样本数 
        }

        return F;
    }
```

总的来说，是利用哈迪温伯格定律来计算的。 1.0 - (hetCount / (2.0 * (1 - altAF) * altAF(double)N ，N是人数。这个值给出的是期望的杂合变异的个数。所以参数F说的就是"实际的hetCount”除以"期望的hetCount"再与1.0取差。当F值越接近0，就意味着实际的hetCount与理论的hetCount越接近。

欢迎关注我的个人公众号：**helixminer（碱基矿工）**

![helixminer-QRCode](https://static.fungenomics.com/images/2021/03/helixminer-mid-red-20210327224600392-20210327224626513.png)
