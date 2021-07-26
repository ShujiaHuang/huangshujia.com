---
title: '三代基因组测序技术原理简介'
date: 2013-08-02 01:00:00+0800
description: 基因组测序技术从第一代到第三代。
image: https://static.fungenomics.com/images/2021/03/02220617-ca30733231724ac4b357b0b98616e898-20210327222526691.png
categories:
    - 生物信息
    - 基因组学
tags:
    - NGS
    - 测序技术

---

## 这篇文章写于2013年，首发在cnblog上，目前已经比较旧了。我重新在WGS系列中对其进行重写，[这里跳转](https://huangshujia.com/2017/08/2017-08-04-begining-wgs-data-analysis-sequecing-tech/)

***

<p><strong>摘要：</strong>从1977年第一代DNA测序技术（Sanger法）<sup>1</sup>，发展至今三十多年时间，测序技术已取得了相当大的发展，从第一代到第三代乃至第四代，测序读长从长到短，再从短到长。虽然就当前形势看来第二代短读长测序技术在全球测序市场上仍然占有着绝对的优势位置，但第三和第四代测序技术也已在这一两年的时间中快速发展着。测序技术的每一次变革，也都对基因组研究，疾病医疗研究，药物研发，育种等领域产生巨大的推动作用。在这里我主要对当前的测序技术以及它们的测序原理做一个简单的小结。</p>

![Figure1](https://static.fungenomics.com/images/2021/03/02214339-89514ce8dd50409392c389f35bb9be01-20210327222547754.png)

<p align="center"><a href="https://images.cnitblog.com/blog/346148/201308/02214339-89514ce8dd50409392c389f35bb9be01.png">图1：测序技术的发展历程</a></p>

<p>生命体遗传信息的快速获得对于生命科学的研究有着十分重要的意义。以上图1（右键打开图片可查看大图，下同）所描述的是自沃森和克里克在1953年建立DNA双螺旋结构以来，整个测序技术的发展历程。</p>
<p><strong><span style="color: #ff0000;">第一代测序技术</span></strong></p>

第一代DNA测序技术用的是1975年由桑格（Sanger）和考尔森（Coulson）开创的链终止法或者是1976-1977年由马克西姆（Maxam）和吉尔伯特（Gilbert）发明的化学法（链降解）. 并在1977年，桑格测定了第一个基因组序列，是噬菌体X174的，全长5375个碱基<sup>1</sup>。自此，人类获得了窥探生命遗传差异本质的能力，并以此为开端步入基因组学时代。研究人员在Sanger法的多年实践之中不断对其进行改进。在2001年，完成的首个人类基因组图谱就是以改进了的Sanger法为其测序基础，Sanger法核心原理是：由于ddNTP的2&rsquo;和3&rsquo;都不含羟基，其在DNA的合成过程中不能形成磷酸二酯键，因此可以用来中断DNA合成反应，在4个DNA合成反应体系中分别加入一定比例带有放射性同位素标记的ddNTP（分为：ddATP,ddCTP,ddGTP和ddTTP），通过凝胶电泳和放射自显影后可以根据电泳带的位置确定待测分子的DNA序列（图2）。这个[网址](https://smcg.cifn.unam.mx/enp-unam/03-EstructuraDelGenoma/animaciones/secuencia.swf)为sanger测序法制作了一个小短片，形象而生动。

<p>值得注意的是，就在测序技术起步发展的这一时期中，除了Sanger法之外还出现了一些其他的测序技术，如焦磷酸测序法、链接酶法等。其中，焦磷酸测序法是后来Roche公司454技术所使用的测序方法<sup>2&ndash;4</sup>，而连接酶测序法是后来ABI公司SOLID技术使用的测序方法<sup>2,4</sup>，但他们的共同核心手段都是利用了Sanger<sup>1</sup>中的可中断DNA合成反应的dNTP。</p>

![Figure2](https://static.fungenomics.com/images/2021/03/02215416-9a87f1e5ab1048a8bbc62cb7605f1bcf-20210327222548001.png)

<p align="center"><a href="https://images.cnitblog.com/blog/346148/201308/02215416-9a87f1e5ab1048a8bbc62cb7605f1bcf.png">图2：Sanger法测序原理</a></p>

<p><strong><span style="color: #ff0000;">第二代测序技术</span></strong></p>
<p>总的说来，第一代测序技术的主要特点是测序读长可达1000bp，准确性高达99.999%，但其测序成本高，通量低等方面的缺点，严重影响了其真正大规模的应用。因而第一代测序技术并不是最理想的测序方法。经过不断的技术开发和改进，以Roche公司的454技术、illumina公司的Solexa，Hiseq技术和ABI公司的Solid技术为标记的第二代测序技术诞生了。第二代测序技术大大降低了测序成本的同时，还大幅提高了测序速度，并且保持了高准确性，以前完成一个人类基因组的测序需要3年时间，而使用二代测序技术则仅仅需要1周，但在序列读长方面比起第一代测序技术则要短很多。表1和图3对第一代和第二代测序技术各自的特点以及测序成本作了一个简单的比较<sup>5</sup>，以下我将对这三种主要的第二代测序技术的主要原理和特点作一个简单的介绍。&nbsp;</p>

![Figure3](https://static.fungenomics.com/images/2021/03/02215522-1d276b88d49a41439f34bb60bf8b7e9f-20210327222548149.png)

<p align="center"><a href="https://images.cnitblog.com/blog/346148/201308/02215522-1d276b88d49a41439f34bb60bf8b7e9f.png">图3. 测序成本的变化</a></p>

<ol>
<li>Illumine</li>
</ol>
<p style="margin-left: 60px;">Illumina公司的Solexa和Hiseq应该说是目前全球使用量最大的第二代测序机器，这两个系列的技术核心原理是相同的<sup>2,4</sup>。这两个系列的机器采用的都是边合成边测序的方法，它的测序过程主要分为以下4步，如图4.</p>
<p>&nbsp;&nbsp;&nbsp;&nbsp; （1）DNA待测文库构建</p>
<p style="margin-left: 60px;">利用超声波把待测的DNA样本打断成小片段，目前除了组装之外和一些其他的特殊要求之外，主要是打断成200-500bp长的序列片段，并在这些小片段的两端添加上不同的接头，构建出单链DNA文库。</p>
<p>&nbsp;&nbsp;&nbsp;&nbsp; （2）Flowcell</p>
<p style="margin-left: 60px;">Flowcell是用于吸附流动DNA片段的槽道，当文库建好后，这些文库中的DNA在通过flowcell的时候会随机附着在flowcell表面的channel上。每个Flowcell有8个channel，每个channel的表面都附有很多接头，这些接头能和建库过程中加在DNA片段两端的接头相互配对（这就是为什么flowcell能吸附建库后的DNA的原因），并能支持DNA在其表面进行桥式PCR的扩增。</p>
<p>&nbsp;&nbsp;&nbsp;&nbsp; （3）桥式PCR扩增与变性</p>
<p style="margin-left: 60px;">桥式PCR以Flowcell表面所固定的接头为模板，进行桥形扩增，如图4.a所示。经过不断的扩增和变性循环，最终每个DNA片段都将在各自的位置上集中成束，每一个束都含有单个DNA模板的很多分拷贝，进行这一过程的目的在于实现将碱基的信号强度放大，以达到测序所需的信号要求。&nbsp;</p>
<p style="margin-left: 30px;">（4）测序</p>
<p style="margin-left: 60px;">测序方法采用边合成边测序的方法。向反应体系中同时添加DNA聚合酶、接头引物和带有碱基特异荧光标记的4中dNTP（如同Sanger测序法）。这些dNTP的3&rsquo;-OH被化学方法所保护，因而每次只能添加一个dNTP。在dNTP被添加到合成链上后，所有未使用的游离dNTP和DNA聚合酶会被洗脱掉。接着，再加入激发荧光所需的缓冲液，用激光激发荧光信号，并有光学设备完成荧光信号的记录，最后利用计算机分析将光学信号转化为测序碱基。这样荧光信号记录完成后，再加入化学试剂淬灭荧光信号并去除dNTP 3&rsquo;-OH保护基团，以便能进行下一轮的测序反应。Illumina的这种测序技术每次只添加一个dNTP的特点能够很好的地解决同聚物长度的准确测量问题，它的主要测序错误来源是碱基的替换，目前它的测序错误率在1%-1.5%之间，测序周期以人类基因组重测序为例，30x测序深度大约为1周。</p>
<p>&nbsp;</p>

![Figure4_a](https://static.fungenomics.com/images/2021/03/02215633-506fdf6d10af4035812ba4898c02985b-20210327222548243.png)
![Figure4_b](https://static.fungenomics.com/images/2021/03/02215656-baeb82e7b0ba4378a7a3417c450243b7-20210327222548523.png)

<p align="center"><a href="https://images.cnitblog.com/blog/346148/201308/02215633-506fdf6d10af4035812ba4898c02985b.png">图4. Illumina测序流程</a></p>

<ol>
<li>Roche 454</li>
</ol>
<p style="margin-left: 60px;">Roche 454测序系统是第一个商业化运营二代测序技术的平台。它的主要测序原理是（图5 abc）<sup>2</sup>：</p>
<p style="margin-left: 30px;">（1）DNA文库制备</p>
<p style="margin-left: 30px;">454测序系统的文件构建方式和illumina的不同，它是利用喷雾法将待测DNA打断成300-800bp长的小片段，并在片段两端加上不同的接头，或将待测DNA变性后用杂交引物进行PCR扩增，连接载体，构建单链DNA文库（图5a）。</p>
<p style="margin-left: 30px;">（2）Emulsion PCR （乳液PCR，其实是一个注水到油的独特过程）</p>
<p style="margin-left: 60px;">454当然DNA扩增过程也和illumina的截然不同，它将这些单链DNA结合在水油包被的直径约28um的磁珠上，并在其上面孵育、退火。</p>
<p style="margin-left: 30px;">乳液PCR最大的特点是可以形成数目庞大的独立反应空间以进行DNA扩增。其关键技术是&ldquo;注水到油&rdquo;（水包油），基本过程是在PCR反应前，将包含PCR所有反应成分的水溶液注入到高速旋转的矿物油表面，水溶液瞬间形成无数个被矿物油包裹的小水滴。这些小水滴就构成了独立的PCR反应空间。理想状态下，每个小水滴只含一个DNA模板和一个磁珠。</p>
<p style="margin-left: 30px;">这些被小水滴包被的磁珠表面含有与接头互补的DNA序列，因此这些单链DNA序列能够特异地结合在磁珠上。同时孵育体系中含有PCR反应试剂，所以保证了每个与磁珠结合的小片段都能独立进行PCR扩增，并且扩增产物仍可以结合到磁珠上。当反应完成后，可以破坏孵育体系并将带有DNA的磁珠富集下来。进过扩增，每个小片段都将被扩增约100万倍，从而达到下一步测序所要求的DNA量。</p>
<p style="margin-left: 30px;">（3）焦磷酸测序</p>
<p style="margin-left: 30px;">测序前需要先用一种聚合酶和单链结合蛋白处理带有DNA的磁珠，接着将磁珠放在一种PTP平板上。这种平板上特制有许多直径约为44um的小孔，每个小孔仅能容纳一个磁珠，通过这种方法来固定每个磁珠的位置，以便检测接下来的测序反应过程。　　</p>
<p style="margin-left: 30px;">测序方法采用焦磷酸测序法，将一种比PTP板上小孔直径更小的磁珠放入小孔中，启动测序反应。测序反应以磁珠上大量扩增出的单链DNA为模板，每次反应加入一种dNTP进行合成反应。如果dNTP能与待测序列配对，则会在合成后释放焦磷酸基团。释放的焦磷酸基团会与反应体系中的ATP硫酸化学酶反应生成ATP。生成的ATP和荧光素酶共同氧化使测序反应中的荧光素分子并发出荧光，同时由PTP板另一侧的CCD照相机记录，最后通过计算机进行光信号处理而获得最终的测序结果。由于每一种dNTP在反应中产生的荧光颜色不同，因此可以根据荧光的颜色来判断被测分子的序列。反应结束后，游离的dNTP会在双磷酸酶的作用下降解ATP，从而导致荧光淬灭，以便使测序反应进入下一个循环。由于454测序技术中，每个测序反应都在PTP板上独立的小孔中进行，因而能大大降低相互间的干扰和测序偏差。454技术最大的优势在于其能获得较长的测序读长，当前454技术的平均读长可达400bp，并且454技术和illumina的Solexa和Hiseq技术不同，它最主要的一个缺点是无法准确测量同聚物的长度，如当序列中存在类似于PolyA的情况时，测序反应会一次加入多个T，而所加入的T的个数只能通过荧光强度推测获得，这就有可能导致结果不准确。也正是由于这一原因，454技术会在测序过程中引入插入和缺失的测序错误。&nbsp;</p>

![Figure5_a](https://static.fungenomics.com/images/2021/03/02215752-3e82bec8a4bd401e928188af86e3f129-20210327222548595.png)
![Figure5_b](https://static.fungenomics.com/images/2021/03/02215812-879851082dfb41e0a2903dbcdedabcd4-20210327222548736.png)
![Figure5_c](https://static.fungenomics.com/images/2021/03/02215833-bfc37705efd74b10be7f1e37e1f851c1-20210327222548812.png)

<p align="center"><a href="https://images.cnitblog.com/blog/346148/201308/02215752-3e82bec8a4bd401e928188af86e3f129.png">图5. Roche 454测序流程</a></p>

<ol>
<li>Solid技术</li>
</ol>
<p style="margin-left: 60px;">Solid测序技术是ABI公司于2007年开始投入用于商业测序应用的仪器。它基于连接酶法，即利用DNA连接酶在连接过程之中测序（图6）<sup>2</sup>,<sup>4</sup>。它的原理是：</p>

![Figure6_a](https://static.fungenomics.com/images/2021/03/02215942-1a81b422e45b40d29f2075656d087391-20210327222548912.png)

<p align="center"><a href="https://images.cnitblog.com/blog/346148/201308/02215942-1a81b422e45b40d29f2075656d087391.png">图6-a. Solid测序技术</a></p>

<p>（1）DNA文库构建</p>
<p>    片段打断并在片段两端加上测序接头，连接载体，构建单链DNA文库。</p>
<p>（2）Emulsion PCR</p>
<p style="margin-left: 60px;">Solid的PCR过程也和454的方法类似，同样采用小水滴emulsion PCR，但这些微珠比起454系统来说则要小得多，只有1um。在扩增的同时对扩增产物的3&rsquo;端进行修饰，这是为下一步的测序过程作的准备。3&rsquo;修饰的微珠会被沉积在一块玻片上。在微珠上样的过程中，沉积小室将每张玻片分成1个、4个或8个测序区域（图6-a）。Solid系统最大的优点就是每张玻片能容纳比454更高密度的微珠，在同一系统中轻松实现更高的通量。</p>
<p>（3）连接酶测序</p>
<p style="margin-left: 60px;">这一步是Solid测序的独特之处。它并没有采用以前测序时所常用的DNA聚合酶，而是采用了连接酶。Solid连接反应的底物是8碱基单链荧光探针混合物，这里将其简单表示为：3&rsquo;-XXnnnzzz-5&rsquo;。连接反应中，这些探针按照碱基互补规则与单链DNA模板链配对。探针的5&rsquo;末端分别标记了CY5、Texas Red、CY3、6-FAM这4种颜色的荧光染料（图6-a）。这个8碱基单链荧光探针中，第1和第2位碱基（XX）上的碱基是确定的，并根据种类的不同在6-8位（zzz）上加上了不同的荧光标记。这是Solid的独特测序法，两个碱基确定一个荧光信号，相当于一次能决定两个碱基。这种测序方法也称之为两碱基测序法。当荧光探针能够与DNA模板链配对而连接上时，就会发出代表第1，2位碱基的荧光信号，图6-a和图6-b中的比色版所表示的是第1，2位碱基的不同组合与荧光颜色的关系。在记录下荧光信号后，通过化学方法在第5和第6位碱基之间进行切割，这样就能移除荧光信号，以便进行下一个位置的测序。不过值得注意的是，通过这种测序方法，每次测序的位置都相差5位。即第一次是第1、2位，第二次是第6、7位&hellip;&hellip;在测到末尾后，要将新合成的链变性，洗脱。接着用引物n-1进行第二轮测序。引物n-1与引物n的区别是，二者在与接头配对的位置上相差一个碱基（图6-a. 8）。也即是，通过引物n-1在引物n的基础上将测序位置往3&rsquo;端移动一个碱基位置，因而就能测定第0、1位和第5、6位&hellip;&hellip;第二轮测序完成，依此类推，直至第五轮测序，最终可以完成所有位置的碱基测序，并且每个位置的碱基均被检测了两次。该技术的读长在2&times;50bp，后续序列拼接同样比较复杂。由于双次检测，这一技术的原始测序准确性高达99.94%，而15x覆盖率时的准确性更是达到了99.999%，应该说是目前第二代测序技术中准确性最高的了。但在荧光解码阶段，鉴于其是双碱基确定一个荧光信号，因而一旦发生错误就容易产生连锁的解码错误。</p>

![Figure6_b](https://static.fungenomics.com/images/2021/03/02220036-548ce076a79149609da6e40075129c4f-20210327222549033.png)

<p align="center"><a href="https://images.cnitblog.com/blog/346148/201308/02220036-548ce076a79149609da6e40075129c4f.png">图6-b. Solid测序技术</a></p>

<p style="margin-left: 30px;"><strong><span style="color: #ff0000;">第三代测序技术</span></strong></p>
<p>测序技术在近两三年中又有新的里程碑。以PacBio公司的SMRT和Oxford Nanopore Technologies纳米孔单分子测序技术，被称之为第三代测序技术。与前两代相比，他们最大的特点就是单分子测序，测序过程无需进行PCR扩增。</p>
<p>其中PacBio SMRT技术其实也应用了边合成边测序的思想<sup>5</sup>，并以SMRT芯片为测序载体。基本原理是： DNA聚合酶和模板结合,4色荧光标记 4 种碱基（即是dNTP）,在碱基配对阶段,不同碱基的加入,会发出不同光,根据光的波长与峰值可判断进入的碱基类型。同时这个 DNA 聚合酶是实现超长读长的关键之一,读长主要跟酶的活性保持有关,它主要受激光对其造成的损伤所影响。PacBio SMRT技术的一个关键是怎样将反应信号与周围游离碱基的强大荧光背景区别出来。他们利用的是ZMW（零模波导孔）原理：如同微波炉壁上可看到的很多密集小孔。小孔直径有考究,如果直径大于微波波长,能量就会在衍射效应的作用下穿透面板而泄露出来，从而与周围小孔相互干扰。如果孔径小于波长,能量不会辐射到周围，而是保持直线状态（光衍射的原理）,从而可起保护作用。同理,在一个反应管(SMRTCell:单分子实时反应孔)中有许多这样的圆形纳米小孔, 即 ZMW(零模波导孔),外径 100多纳米,比检测激光波长小(数百纳米),激光从底部打上去后不能穿透小孔进入上方溶液区,能量被限制在一个小范围(体积20X 10-21 L)里,正好足够覆盖需要检测的部分,使得信号仅来自这个小反应区域,孔外过多游离核苷酸单体依然留在黑暗中,从而实现将背景降到最低。另外，可以通过检测相邻两个碱基之间的测序时间，来检测一些碱基修饰情况，既如果碱基存在修饰，则通过聚合酶时的速度会减慢，相邻两峰之间的距离增大，可以通过这个来之间检测甲基化等信息（图7）。SMRT技术的测序速度很快，每秒约10个dNTP。但是，同时其测序错误率比较高（这几乎是目前单分子测序技术的通病），达到15%,但好在它的出错是随机的，并不会像第二代测序技术那样存在测序错误的偏向，因而可以通过多次测序来进行有效的纠错。</p>

![Figure7](https://static.fungenomics.com/images/2021/03/02220111-aac12c31ed944097b1ef793332ed5b24-20210327222549177.png)

<p align="center"><a href="https://images.cnitblog.com/blog/346148/201308/02220111-aac12c31ed944097b1ef793332ed5b24.png">图7. PacBio SMRT测序原理</a></p>

<p>Oxford Nanopore Technologies公司所开发的纳米单分子测序技术与以往的测序技术皆不同，它是基于电信号而不是光信号的测序技术<sup>5</sup>。该技术的关键之一是，他们设计了一种特殊的纳米孔，孔内共价结合有分子接头。当DNA碱基通过纳米孔时，它们使电荷发生变化，从而短暂地影响流过纳米孔的电流强度（每种碱基所影响的电流变化幅度是不同的），灵敏的电子设备检测到这些变化从而鉴定所通过的碱基（图8）。</p>
<p>该公司在去年基因组生物学技术进展年会(AGBT)上推出第一款商业化的纳米孔测序仪，引起了科学界的极大关注。纳米孔测序（和其他第三代测序技术）有望解决目前测序平台的不足，纳米孔测序的主要特点是：读长很长，大约在几十kb，甚至100 kb;错误率目前介于1%至4%，且是随机错误，而不是聚集在读取的两端;数据可实时读取;通量很高(30x人类基因组有望在一天内完成);起始DNA在测序过程中不被破坏;以及样品制备简单又便宜。理论上，它也能直接测序RNA。</p>
<p>纳米孔单分子测序计算还有另一大特点，它能够直接读取出甲基化的胞嘧啶，而不必像传统方法那样对基因组进行bisulfite处理。这对于在基因组水平直接研究表观遗传相关现象有极大的帮助。并且改方法的测序准确性可达99.8%，而且一旦发现测序错误也能较容易地进行纠正。但目前似乎还没有应用该技术的相关报道。</p>

![Figure8](https://static.fungenomics.com/images/2021/03/02220146-57c9245f4a1b4f0bbd3a0caf7046f7bc-20210327222549319.png)

<p align="center"><a href="https://images.cnitblog.com/blog/346148/201308/02220146-57c9245f4a1b4f0bbd3a0caf7046f7bc.png">图8. 纳米孔测序</a></p>

<p>其他测序技术</p>
<p>目前还有一种基于半导体芯片的新一代革命性测序技术&mdash;&mdash;Ion Torrent<sup>6</sup>。该技术使用了一种布满小孔的高密度半导体芯片， 一个小孔就是一个测序反应池。当DNA聚合酶把核苷酸聚合到延伸中的DNA链上时，会释放出一个氢离子，反应池中的PH发生改变，位于池下的离子感受器感受到H+离子信号，H+离子信号再直接转化为数字信号，从而读出DNA序列（图9）。这一技术的发明人同时也是454测序技术的发明人之一&mdash;&mdash;Jonathan Rothberg，它的文库和样本制备跟454技术很像，甚至可以说就是454的翻版，只是测序过程中不是通过检测焦磷酸荧光显色，而是通过检测H+信号的变化来获得序列碱基信息。Ion Torrent相比于其他测序技术来说，不需要昂贵的物理成像等设备，因此，成本相对来说会低，体积也会比较小，同时操作也要更为简单，速度也相当快速，除了2天文库制作时间，整个上机测序可在2-3.5小时内完成，不过整个芯片的通量并不高，目前是10G左右，但非常适合小基因组和外显子验证的测序。&nbsp;&nbsp;&nbsp;&nbsp;</p>

![Figure9_a](https://static.fungenomics.com/images/2021/03/02220244-63a58f22939d4b59890f8c95068a4296-20210327222549400.png)
![Figure9_b](https://static.fungenomics.com/images/2021/03/02220317-7bd9cecbe1c04cdbad13a63e66eeafe5-20210327222549472.png)
![Figure9_c](https://static.fungenomics.com/images/2021/03/02220354-9cd0104936984b47b19a64f0b8ed0d68-20210327222549719.png)

<p align="center"><a href="https://images.cnitblog.com/blog/346148/201308/02220244-63a58f22939d4b59890f8c95068a4296.png">图9. Ion Torrent</a></p>

<p align="center">&nbsp;</p>
<p><strong>小结</strong></p>
<p>以上，对各代测序技术的原理做了简要的阐述，这三代测序技术的特点比较汇总在以下表1和表2中。其中测序成本，读长和通量是评估该测序技术先进与否的三个重要指标。第一代和第二代测序技术除了通量和成本上的差异之外，其测序核心原理（除Solid是边连接边测序之外）都是基于边合成边测序的思想。第二代测序技术的优点是成本较之一代大大下降，通量大大提升，但缺点是所引入PCR过程会在一定程度上增加测序的错误率，并且具有系统偏向性，同时读长也比较短。第三代测序技术是为了解决第二代所存在的缺点而开发的，它的根本特点是单分子测序，不需要任何PCR的过程，这是为了能有效避免因PCR偏向性而导致的系统错误，同时提高读长，并要保持二代技术的高通量，低成本的优点。</p>

<p align="left">表1：测序技术的比较</p>
<div align="center">
<table style="width: 727px; height: 1267px;" border="1" cellspacing="0" cellpadding="0">
<tbody>
<tr>
<td width="48">
<p align="center"><strong>第</strong><strong>X</strong><strong>代</strong></p>
</td>
<td width="48">
<p align="center"><strong>公司</strong></p>
</td>
<td width="64">
<p align="center"><strong>平台名称</strong></p>
</td>
<td width="59">
<p align="center"><strong>测序方法</strong></p>
</td>
<td width="53">
<p align="center"><strong>检测方法</strong></p>
</td>
<td width="56">
<p align="center"><strong>大约读长</strong><strong>(</strong><strong>碱基数</strong><strong>)</strong></p>
</td>
<td width="91">
<p align="center"><strong>优点</strong></p>
</td>
<td width="112">
<p align="center"><strong>相对局限性</strong></p>
</td>
</tr>
<tr>
<td width="48">
<p align="center">第一代</p>
</td>
<td width="48">
<p align="center">ABI/生命技术公司</p>
</td>
<td width="64">
<p align="center">3130xL-3730xL</p>
</td>
<td width="59">
<p align="center">桑格-毛细管电泳测序法</p>
</td>
<td width="53">
<p align="center">荧光/光学</p>
</td>
<td width="56">
<p align="center">600-1000</p>
</td>
<td width="91">
<p align="left">高读长，准确度一次性达标率高，能很好处理重复序列和多聚序列</p>
</td>
<td width="112">
<p align="left">通量低；样品制备成本高，使之难以做大量的平行测序</p>
</td>
</tr>
<tr>
<td width="48">
<p align="center">第一代</p>
</td>
<td width="48">
<p align="center">贝克曼</p>
</td>
<td width="64">
<p align="center">GeXP遗传分析系统</p>
</td>
<td width="59">
<p align="center">桑格-毛细管电泳测序法</p>
</td>
<td width="53">
<p align="center">荧光/光学</p>
</td>
<td width="56">
<p align="center">600-1000</p>
</td>
<td width="91">
<p align="left">高读长，准确度一次性达标率高，能很好处理重复序列和多聚序列；易小型化</p>
</td>
<td width="112">
<p align="left">通量低；单个样品的制备成本相对较高</p>
</td>
</tr>
<tr>
<td width="48">
<p align="center">第二代</p>
</td>
<td width="48">
<p align="center">Roche/454</p>
</td>
<td width="64">
<p align="center">基因组测序仪FLX系统</p>
</td>
<td width="59">
<p align="center">焦磷酸测序法</p>
</td>
<td width="53">
<p align="center">光学</p>
</td>
<td width="56">
<p align="center">230-400</p>
</td>
<td width="91">
<p align="left">在第二代中最高读长；比第一代的测序通量大</p>
</td>
<td width="112">
<p align="left">样品制备较难；难于处理重复和同种碱基多聚区域；试剂冲洗带来错误累积；仪器昂贵</p>
</td>
</tr>
<tr>
<td width="48">
<p align="center">第二代</p>
</td>
<td width="48">
<p align="center">Illumina</p>
</td>
<td width="64">
<p align="center">HiSeq2000,HiSeq2500/MiSeq</p>
</td>
<td width="59">
<p align="center">可逆链终止物和合成测序法</p>
</td>
<td width="53">
<p align="center">荧光/光学</p>
</td>
<td width="56">
<p align="center"><strong>2x150</strong></p>
</td>
<td width="91">
<p align="left">很高测序通量</p>
</td>
<td width="112">
<p align="left">仪器昂贵；用于数据删节和分析的费用很高</p>
</td>
</tr>
<tr>
<td width="48">
<p align="center">第二代</p>
</td>
<td width="48">
<p align="center">ABI/Solid</p>
</td>
<td width="64">
<p align="center">5500xlSolid系统</p>
</td>
<td width="59">
<p align="center">连接测序法</p>
</td>
<td width="53">
<p align="center">荧光/光学</p>
</td>
<td width="56">
<p align="center">25-35</p>
</td>
<td width="91">
<p align="left">很高测序通量；在广为接受的几种第二代平台中，所要拼接出人类基因组的试剂成本最低</p>
</td>
<td width="112">
<p align="left">测序运行时间长；读长短，造成成本高，数据分析困难和基因组拼接困难；仪器昂贵</p>
</td>
</tr>
<tr>
<td width="48">
<p align="center">第二代</p>
</td>
<td width="48">
<p align="center">赫利克斯</p>
</td>
<td width="64">
<p align="center">Heliscope</p>
</td>
<td width="59">
<p align="center">单分子合成测序法</p>
</td>
<td width="53">
<p align="center">荧光/光学</p>
</td>
<td width="56">
<p align="center">25-30</p>
</td>
<td width="91">
<p align="left">高通量；在第二代中属于单分子性质的测序技术</p>
</td>
<td width="112">
<p align="left">读长短，推高了测序成本，降低了基因组拼接的质量；仪器非常昂贵</p>
</td>
</tr>
<tr>
<td width="48">
<p align="center">第三代</p>
</td>
<td width="48">
<p align="center">太平洋生物科学公司</p>
</td>
<td width="64">
<p align="center">PacBio RS</p>
</td>
<td width="59">
<p align="center">实时单分子DNA测序</p>
</td>
<td width="53">
<p align="center">荧光/光学</p>
</td>
<td width="56">
<p align="center">~1000</p>
</td>
<td width="91">
<p align="left">高平均读长，比第一代的测序时间降低；不需要扩增；最长单个读长接近3000碱基</p>
</td>
<td width="112">
<p align="left">并不能高效地将DNA聚合酶加到测序阵列中；准确性一次性达标的机会低（81-83%）；DNA聚合酶在阵列中降解；总体上每个碱基测序成本高（仪器昂贵）；</p>
</td>
</tr>
<tr>
<td width="48">
<p align="center">第三代</p>
</td>
<td width="48">
<p align="center">全基因组学公司</p>
</td>
<td width="64">
<p align="center">GeXP遗传分析系统</p>
</td>
<td width="59">
<p align="center">复合探针锚杂交和连接技术</p>
</td>
<td width="53">
<p align="center">荧光/光学</p>
</td>
<td width="56">
<p align="center">10</p>
</td>
<td width="91">
<p align="left">在第三代中通量最高；在所有测序技术中，用于拼接一个人基因组的试剂成本最低；每个测序步骤独立，使错误的累积变得最低</p>
</td>
<td width="112">
<p align="left">低读长；&nbsp;模板制备妨碍长重复序列区域测序；样品制备费事；尚无商业化供应的仪器</p>
</td>
</tr>
<tr>
<td width="48">
<p align="center">第三代</p>
</td>
<td width="48">
<p align="center">Ion Torrent/生命技术公司</p>
</td>
<td width="64">
<p align="center">个人基因组测序仪（PGM）</p>
</td>
<td width="59">
<p align="center">&nbsp;合成测序法</p>
</td>
<td width="53">
<p align="center">以离子敏感场效应晶体管检测pH值变化</p>
</td>
<td width="56">
<p align="center">100-200</p>
</td>
<td width="91">
<p align="left">对核酸碱基的掺入可直接测定；在自然条件下进行DNA合成（不需要使用修饰过的碱基）</p>
</td>
<td width="112">
<p align="left">一步步的洗脱过程可导致错误累积；阅读高重复和同种多聚序列时有潜在困难；</p>
</td>
</tr>
<tr>
<td width="48">
<p align="center">第三代</p>
</td>
<td width="48">
<p align="center">牛津纳米孔公司</p>
</td>
<td width="64">
<p align="center">&nbsp;gridION</p>
</td>
<td width="59">
<p align="center">纳米孔外切酶测序</p>
</td>
<td width="53">
<p align="center">电流</p>
</td>
<td width="56">
<p align="center">尚未定量</p>
</td>
<td width="91">
<p align="left">有潜力达到高读长；可以成本生产纳米孔；无需荧光标记或光学手段</p>
</td>
<td width="112">
<p align="left">切断的核苷酸可能被读错方向；难于生产出带多重平行孔的装置</p>
</td>
</tr>
</tbody>
</table>
</div>
<p align="center">&nbsp;&nbsp;</p>
<p align="left">表2：主流测序机器的成本测序比较</p>

![Table](https://static.fungenomics.com/images/2021/03/02220513-79289588a8bf40e081f3ab8c9cdda24d-20210327222549896.png)

以下图10展示了当前全球测序仪的分布情况。图中的几个热点区主要分布在中国的深圳（主要是华大），南欧，西欧和美国。

![Figure10](https://static.fungenomics.com/images/2021/03/02220617-ca30733231724ac4b357b0b98616e898-20210327222550124.png)

<p align="center"><a href="https://omicsmaps.com/">图10. 测序仪全球分布</a></p>

*参考文献*

<blockquote>
<p>1. Sanger, F. &amp; Nicklen, S. DNA sequencing with chain-terminating. <strong>74</strong>, 5463&ndash;5467 (1977).</p>
<p>2. Mardis, E. R. Next-generation DNA sequencing methods. <em>Annual review of genomics and human genetics</em> <strong>9</strong>, 387&ndash;402 (2008).</p>
<p>3. Shendure, J. &amp; Ji, H. Next-generation DNA sequencing. <em>Nature biotechnology</em> <strong>26</strong>, 1135&ndash;45 (2008).</p>
<p>4. Metzker, M. L. Sequencing technologies - the next generation. <em>Nature reviews. Genetics</em> <strong>11</strong>, 31&ndash;46 (2010).</p>
<p>5. Niedringhaus, T. P., Milanova, D., Kerby, M. B., Snyder, M. P. &amp; Barron, A. E. Landscape of Next-Generation Sequencing Technologies. 4327&ndash;4341 (2011).</p>
<p>6. Rothberg, J. M. <em>et al.</em> An integrated semiconductor device enabling non-optical genome sequencing. <em>Nature</em> <strong>475</strong>, 348&ndash;52 (2011).&nbsp;</p>
</blockquote>




--------------

欢迎关注我的个人公众号：**helixminer（碱基矿工）**
![helixminer-QRCode](https://static.fungenomics.com/images/2021/03/helixminer-mid-red-20210327222459389-20210327222550461.png)