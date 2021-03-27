---
title: '目前最好最完整的SOAPdenovo使用说明'
date: 2015-07-09 
image: http://image.fungenomics.com/tree.jpeg
description: 从0到1的基因组序列组装总是最难的。
categories:
    - 生物信息
tags:
    - 组装
    - SOAPdenovo
---

由于丹麦人国家基因组项目的原因，近期我整理了一份关于SOAPdenovo2的使用说明，内容包括了程序使用、参数的详细说明、参数如何调整、各个主要输出文件的格式说明等，而且我敢说这是目前最好最全的！

### 简介
SOAPdenovo（目前最新版是SOAPdenovo2）是一种应用de Bruijn graph组装短read的方法，它以kerm为节点单位，利用de Bruijn图的方法实现全基因组的组装，与其他短序列组装软件相比，它可以进行大型基因组，比如人类基因组的组装，组装结果更加准确可靠，可以通过组装的结果非常准确地鉴别出基因组上的序列结构性变异，为构建全基因组参考序列和以低测序成本对未知基因组实施精确分析创造了可能。

下载地址：<http://soap.genomics.org.cn/soapdenovo.html>

安装：
 * 下载SOAPdenovo的压缩包          
 * 解压缩     
 * 将得到可执行文件SOAPdenovo和一个配置文件的模板example.contig

### 使用程序及参数

SOAPdenovo可以一步跑完，也可以分成四步单独跑，一步跑完的脚本:

```bash
./SOAPdenovo all -s lib.cfg -K 29 -D 1 -o ant >>ass.log
```

四步单独跑的脚本:
```bash
./SOAPdenovo pregraph -s lib.cfg -d 1  -K 29 -o ant >pregraph.log
./SOAPdenovo contig -g ant -D 1 -M 3 >contig.log
./SOAPdenovo map -s lib23.cfg -g ant >map.log
./SOAPdenovo scaff -g ant -F >scaff.log
```

### 参数说明

```
用法：/PathToProgram/SOAPdenovo all -s configFile [-K kmer -d KmerFreqCutOff -D EdgeCovCutoff -M mergeLevel -R -u -G gapLenDiff -L minContigLen -p n_cpu] -o Output
```

```
    -s    STR     配置文件
    -o    STR     输出文件的文件名前缀
    -g    STR     输入文件的文件名前缀
    -K    INT     输入的K-mer值大小，默认值23，取值范围 13-63
    -p    INT     程序运行时设定的线程数，默认值8
    -R            利用read鉴别短的重复序列，默认值不进行此操作
    -d    INT     去除频数不大于该值的k-mer，默认值为0
    -D    INT     去除频数不大于该值的由k-mer连接的边，默认值为1，即该边上每个点的频数都小于等于1时才去除
    -M    INT     连接contig时合并相似序列的等级，默认值为1，最大值3。
    -F            利用read对scaffold中的gap进行填补，默认不执行
    -u            构建scaffold前不屏蔽高覆盖度的contig，这里高频率覆盖度指平均contig覆盖深度的2倍。默认屏蔽
    -G    INT     估计gap的大小和实际补gap的大小的差异，默认值为50bp。
    -L            用于构建scaffold的contig的最短长度，默认为：Kmer参数值 ×2
```

### 使用方法及示例

（1）示例
```bash
SOAPdenovo all -s HCB.lib -K 25 -d -o test
```

（2） 输入文件
configFile，配置文件内容如下，非程序生成，需要软件使用者自己配置。各个说明参考如下：
```
# 以“#”开头的行是注释内容

# maximal read length （read的最大长度）
# 该值一般设置的比实际read读长稍微短一些，截去测序最后的部分，具体长度看测序质量
max_rd_len=50  

[LIB] # 文库信息以此开头
# 文库平均插入长度，一般取插入片段分布图中给出的文库大小
avg_ins=200

#序列是否需要被反转，目前的测序技术，插入片段大于等于2k的采用了环化，所以对于插入长度大于等于2k文库，序列需要反转，reverse_seq＝1，小片段设为0
reverse_seq=0

# 该文库中的read序列在组装的哪些过程（contig/scaff/fill）中用到
# 设为1：只用于构建contig；
# 设为2：只用于构建scaffold；
# 设为3：同时用于构建contig和scaffold；
# 设为4：只用于补洞

# [注意]短插入片段(<2K)的设为3，同时用于构建contig和scaffold，长插入片段(>=2k)设为2，不用于构建contig，只用于构建scaffold，454single 长reads只用于补洞。
asm_flags=3

# rank该值取整数，决定了reads用于构建scaffold的次序，值越低，数据越优先用于构建scaffold。
# 设置了同样rank的文库数据会同时用于组装scaffold。
# 一般将短插入片段设为1；2k设为2；5k设为3；10k设为4；
# 当某个档的数据量较大时，也可以将其分为多个档，同样，当某档数据量不足够时，可以将多个档的数据合在一起构建scaffold。
# 这里说的数据量够与不够是从该档的测序覆盖度和物理覆盖度两个方面来考虑的。
rank=1

# 可选参数，pair_num_cutoff该参数规定了连接两个contig 或者是pre-scaffold 的可信连接的阈值，即，当连接数大于该值，连接才算有效。短插入片段(<2k)默认值为3，长插入长度序列默认值为5
pair_num_cutoff=3

# map_len该参数规定了在map过程中 reads和contig的比对长度必须达到该值（比对不容mismacth和gap），该比对才能作为一个可信的比对。可选参数，短插入片段(<2k)一般设置为32，长插入片段设置为35，默认值是K＋2.
map_len=32

# read 1的fastq格式的序列文件，“/path/**LIBNAMEA**/fastq_read_1.fq”为read的存储路径
q1=/path/**LIBNAMEA**/fastq_read_1.fq

# read 2的fastq格式的序列文件，与read1对应的read2文件紧接在read1之后）
q2=/path/**LIBNAMEA**/fastq_read_2.fq

#read 1的fasta格式的序列文件
f1=/path/**LIBNAMEA**/fasta_read_1.fa

# read 2的fasta格式的序列文件
f2=/path/**LIBNAMEA**/fasta_read_2.fa

# 单向测序得到的fastq格式的序列文件
q=/path/**LIBNAMEA**/fastq_read_single.fq

# 单向测序得到的fasta格式的序列文件
f=/path/**LIBNAMEA**/fasta_read_single.fa

# 双向测序得到的一个fasta格式的序列文件
p=/path/**LIBNAMEA**/pairs_in_one_file.fa

```

### 输出文件及说明

SOAPdenovo 分四部分别对应的输出文件：
```bash
 1. pregraph  生成7个文件 *.kmerFreq  *.edge  *.preArc  *.markOnEdge  *.path *.vertex  *.preGraphBasic
 2. contig       生成4个文件 *.contig  *.ContigIndex  *.updated.edge  *.Arc
 3. map          生成3个文件 *.readOnContig  *.peGrads  *.readInGap
 4. scaff        生成6个文件 *.newContigIndex  *.links  *.scaf  *.scaf_gap  *.scafSeq  *.gapSeq
```

*.contig：contig序列文件，fasta格式；

*.scafSeq：fasta格式的scaffold序列文件，contig之间的gap用N填充；

对于得到的*.scafSeq文件还需要用GapCloser去合并其中的gap，最后的contig文件则是对补洞之后的scaffold文件通过打断N区的方法得到。

以上两个文件是组装结果中最主要的输出。

*.scaf：包括scaffold中contig的详细信息；在scaffold行中包括scaffold名字、contig长度和该scaffold长度。在contig行包括contig名字、contig在scaffold上的起始位置、正反链、长度和contig间的链接信息;

*.links：contig间的pair-end连接信息;

*.readOnContig：reads在contig上的位置;

*.peGrads： 主要可以通过调整本文件中的参数来显示构建scaffold所用到的插入片段库的个数，总共要到的read数，最长的read的长度，每个库对应的哪些reads，rank设置，pair_num_cutoff设置。例如：

```bash
grads&num: 10   522083934       70
323     104577616       1       3
334     180770522       1       3
345     226070520       1       3
486     361955834       2       3
2200    392088076       3       5
2290    422272580       3       5
2400    445522690       3       5
4870    475666064       4       5
9000    511030930       5       8
9110    522083934       5       5
```

该文件中共分成4列。组装的配置文件中有n个文库，该文件则有n+1行，且按照文库大小顺序排列。
第1行中，第二三四列分别是 所用文库，reads总数和组装中用到的最长的reads长度。
第2行中，四列分别是文库大小，文库中的reads数目，该文库reads用到的rank等级和该文库中reads用到的pair_num_cutoff。
第3～n+1行，四列分别是文库大小，文库中的reads数目加上前面的文库中的reads总数，该文库reads用到的rank等级和该文库中reads用到的pair_num_cutoff。
如果配置文件中没有设置pair_num_cutoff，即使用默认参数，则最后一列显示为0。

对于SOAPdenovo的每个步骤都有日志文件输出，要保存好日志文件，日志文件中包含有很多有用的信息。

### SOAPdenovo日志输出说明

1）pregraph.log:  其中有很多的统计信息，包括构建debruijn-graph时用到多少reads数，构图中生成了多少uniq的kmer以及设置-d参数后去除了多少kmer。
在pregraph中，可选参数有 –R –K –d  结果如：

```bash
5467781332 nodes allocated, 70662750348 kmer in reads, 70662750348 kmer processed
3283081670 kmer removed
```

其中Kmer 数是取决于所设k值大小以及数据量，nodes数即特异性的kmer数目，当nodes数目过高（一般和基因组大小差不多大小），可能是数据的错误率比较高，也可能是存在杂合。若nodes数目偏小，并且kmer数目很多，则基因组本身可能存在一定的重复度。对于k值的选取，当数据量充足时（>=40X），植物基因组一般采用大kmer会有比较好的效果，而对于动物基因组，k值一般多取27和29则足够。kmer removed表示的 –d 参数所去除的低频的kmer。

2）contig.log: contig 中，可选参数 –R –D –M，注，-R 参数的选定，必须pregraph和contig中同时选择才有效。结果例子：

```
16430183 pairs found, 2334584 pairs of paths compared, 1674493 pairs merged
```

从merged的数量可以作为估计杂合以及测序错误的程度。

```bash
sum up 1932549703bp, with average length 1170
the longest is 36165bp, contig N50 is 2871 bp,contig N90 is 553 bp
```

3）map.log: 
```
Output 415219610 out of 1956217742 (21.2)% reads in gaps
1661094582 out of 1956217742 (84.9)% reads mapped to contigs
```

一般情况下，reads in gap的比例和map to contig 的比例总和大于1。可能是因为reads map到多个地方都被算在其中的原因。当map to contig的比例很高（80%左右时），但是组装效果并不很好，可能是重复序列比较多。reads in gap比例较高（大于40%），是因为基因组组装的较碎，gap区域较多。
map_len 默认值=K+5，当默认值大于设置的map_len时，以默认值为准，当默认值小于map_len值时，设置的map_len为准。

4）scaff.log:
```
average contig coverage is 23, 5832270 contig masked
```

构建scaffold是对高频覆盖的contig进行屏蔽（即频率高于average contig coverage的两倍的contig不用于构建scaffold），从这里可以看出组装的基因组一定的重复情况。

```
estimated PE size 162, by 40034765 pairs
on contigs longer than 173, 38257479 pairs found,SD=8, insert_size estimated: 163
```

173 是配置文件中该文库的insertsize，163 是根据reads  map到contig上的距离的估计值，8是这个分布的标准偏差。一般考虑 比对上去的pair数目和SD值。若pair对数很多且SD值很小（小片段文库数据不超过三位数，大片段文库数据部超过500），那我们一般可以将配置文件中的文库插入片段的值改对短插入片段文库（<1k）的大小估计值，一般是比较准确的，下次组装以及补洞时应根据这个值对原来配置文件中的insertsize信息做修正。对于大片段文库（>=2K），因为是把reads map到contig上，若最长contig较短时，可能找不到成pair比对上去的reads，这时，无法估计文库大小，需要自己将大片段一级一级的map到前一级的组装结果上，然后再分析大片段文库的插入片段大小。注，需要调整insertsize信息时，只需要修改* .peGrads文件中的第一列，然后删除*.links文件，重新跑scaff这一步即可。即构建scaffold时，主要是根据*.links文件的信息进行连接。

```
Cutoff for number of pairs to make a reliable connection: 3
1124104 weak connects removed (there were 4773564 active cnnects))
```

Cutoff for number是在配置文件中设的pair_num_cutoff值，weak connects是低于这个值被认定为无效的连接数，active connects是满足cutoff的连接数，根据这个数值可对pair_num_cutoff做调整

```
Picked  25241 subgraphs,4 have conflicting connections
```

conflicting connections 是表示构建scaffold时的矛盾数，矛盾数比较高（>100）时，可根据前面的有效连接数，适当提高pair_num_cutoff值，即提高scaffold连接要求的最少关系数

```
182483 scaffolds&singleton sum up 1990259817bp, with average length 10906
the longest is 6561520bp,scaffold N50 is 836795 bp, scaffold N90 is 157667 bp
```

scaffold 统计信息，将是根据rank分梯度的统计:
```
Done with 13301 scaffolds, 2161915 gaps finished, 2527441 gaps overall
```

-F 参数补洞的统计信息。

### 参数调整

一般组装时需要调整的参数，主要分两种：

一种是针对脚本中的参数改动：如调整  -K  -R  -d  -D  -M
-K 值一般与基因组的特性和数据量相关，目前用到的SOAPdenovo软件主要有两个版本，grape1123和grape63mer，其中grape1123是最新版的组装软件，K值范围13-31，grape63mer是可以使用大kmer的组装版本，K值范围13-63。 

【经验】：植物基因组的组装采用大kmer效果会比较好（要求短片段reads长度75bp），动物基因组很少有用到大kmer后有明显改进效果的，且动物基因组的组装K值一般设置为27和29较多。

-R参数，对于动物基因组，R参数一般不设置，植物基因组由于较多的repeat区，则设置R参数后，效果更好。注意，设置-R时，一般使用-M 的默认值。（熊猫基因组组装时得出的结论）

-M 参数，0-3,默认值1。一般杂合率为千分之几就设为几。熊猫基因组组装时-M 2 。

-d 参数，对于没有纠错，没有处理的质量又较差的原始数据，kmer的频数为1的很多的数据的组装，一般设置为-d 1 则足够。对于处理过，或者是测序质量较好的数据，可以不用设置。数据量很多时，也可以以-d 参数去除部分质量稍差的数据。

-D 参数，默认为1，一般不用另行设置。

第二种，从map这一过程去调节参数。可以调整配置文件的map_len的值和调整文件*.peGrads。

当文库插入片段分布图中文库大小与实验给出的文库大小差异很大时，调整*.peGrads文件中的插入片段大小。

根据每一档数据的数据量去调整文库的rank等级。当该文库的数据量很多或者是在构建scaffold的过程中的冲突数很多时，可是适当的调大第四列 的pair_num_cutoff，把条件设置的更严一些。

### 内存估计

SOAPdenovo的四个步骤消耗的内存是不一样的，其中第一步消耗的内存最多，使用没有纠错的的reads，(K<=31)第一步消耗的内存在基因组大小的80－100倍左右，纠过错则在40－50倍左右，第二步相对消耗的内存会少很多，第三步消耗的内存是仅次于第一步的，在第一步的一半左右，第四步消耗的内存也会比较少。对于CPU的使用，默认是8个，如果申请内存时申请一个计算节点的所有内存测将CPU就设置为该计算节点的CPU个数充分利用计算资源，如果仅申请一个节点的部分内存则根据实际情况考虑。对于大kemr(K>31)其内存使用是(k<=31)的1.5倍左右，有时甚至更多，要充分估计内存的使用，在第一次运行的时候考虑不能太保守。

### 常见错误

1）配置文件中read存储路径错误

只输出日志文件。
pregraph.log中的错误信息：“Cannot open /path/**LIBNAMEA**/fastq_read_1.fq. Now exit to system...”

2）-g 后所跟参数与pregraph（第一步） -o  后所跟参数名不一致

contig  map  scaff 这三个步骤都只是输出日志文件。

contig.log中的错误信息：“Cannot open *.preGraphBasic. Now exit to system...”

map.log中的错误信息：“Cannot open *.contig. Now exit to system...”

scaff.log中的错误信息：“Cannot open *.preGraphBasic. Now exit to system...”

3）从map开始重新跑时，需要删除*.links文件，否则会生成core文件，程序退出。


欢迎关注我的个人公众号：**helixminer（碱基矿工）**

![helixminer-QRCode](https://static.fungenomics.com/images/2021/03/helixminer-mid-red.png)



