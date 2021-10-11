# famCircle 基因家族分析，tandem可视化分析

## 分析流程
![分析流程](https://images.gitee.com/uploads/images/2021/1011/164430_44551292_8074509.png "famCircle.vpd.png")
>基因组下载以及处理  
>路线：  
>>前期处理序列比对、共线性扫描、ks计算  
>>hmmer 隐马尔可夫模型搜索筛选  
>>screen tandem数量分部  
>>>断开处：  
>>>>pfam结构域筛选# http://pfam.xfam.org/  
>>>>α/β筛选# http://cb.csail.mit.edu/cb/paircoil2/  
>>主程序选择：  
>>>inner  
>>>>typing 数据格式化  
>>>>inner绘图  
>>>outer  
>>>>typing 数据格式化  
>>>>outer绘图  
>>part局部tandem绘制  
>其他可视化分析工具：  
>>ks展示  
>>>全基因组 ks展示  
>>>block ks展示  
>>>circle全基因组共线性展示  
>>>line局部共线性绘制  

# 软件安装与使用
## 下载安装
### 从PyPI:https://pypi.org/project/famCircle/ 下载安装包到本地安装
&ensp;&ensp;&ensp;&ensp;命令  
```
pip install xx/xx/famCircle-*.*.*.tar.gz
pip install xx/xx/famCircle-*.*.*-py3-none-any.whl
```
### 直接通过pip安装
&ensp;&ensp;&ensp;&ensp;命令 

```
pip install famCircle
```
## 软件配置
安装路径下/famCircle/conf.ini

```
[ini]
hmmer_path = ~/bin/
# hmm安装路径
clustalw_path = ~/bin/clustalw
#clustalw安装路径
muscle_path = ~/bin/
#muscle安装路径
mafft_path = ~/bin/
#mafft安装路径
yn00_path = ~/bin/yn00
#paml安装路径
pal2nal_path = ~/bin/pal2nal.pl
#pal2nal安装路径
```
## 软件使用
交流群![交流群号](https://images.gitee.com/uploads/images/2021/0928/123500_6eb999c4_8074509.jpeg "B5726647-A6DC-4C3E-ADF7-30720901BC51.jpeg")
### help

```
famCircle -h
```
### 文件配置 以outer为例

```
famCircle -o ? > total.conf# 覆盖
famCircle -o ? >> total.conf# 追加
```

### 功能使用 以outer为例

```
famCircle -o total.conf
```

# 配置文件解释

## 前期
### 计算ks
&ensp;&ensp;&ensp;&ensp;[Ks]  
&ensp;&ensp;&ensp;&ensp;cds_file = 	cds file# 基因组cds文件  
&ensp;&ensp;&ensp;&ensp;pep_file = 	pep file# 基因组蛋白质文件(可选)  
&ensp;&ensp;&ensp;&ensp;align_software = muscle# 比对软件  
&ensp;&ensp;&ensp;&ensp;pairs_file = gene pairs file# 共线性文件  
&ensp;&ensp;&ensp;&ensp;ks_file = ks result# 输出文件  


### 基因块ks数据分布图
&ensp;&ensp;&ensp;&ensp;[Ks_block]  
&ensp;&ensp;&ensp;&ensp;species_list = specise name# 物种列表例如：ath,vvi  
&ensp;&ensp;&ensp;&ensp;ks = ks file# ks文件  
&ensp;&ensp;&ensp;&ensp;area = 0,2.5# ks绘制区间  
&ensp;&ensp;&ensp;&ensp;model = YN00/NG86# 绘制ks来源  
&ensp;&ensp;&ensp;&ensp;blockfile = block file# 共线性文件  
&ensp;&ensp;&ensp;&ensp;blocklength = 6# block长度控制参数  
&ensp;&ensp;&ensp;&ensp;pvalue = 1e-5# block阈值控制参数  
&ensp;&ensp;&ensp;&ensp;savecsv = save csv# 筛选结果文件  
&ensp;&ensp;&ensp;&ensp;savefile = save file (*.png, *.pdf)# 筛选结果可视化  
例：![水稻block的ks分布](https://images.gitee.com/uploads/images/2021/1011/164725_6a93e1f7_8074509.png "osa_osa.new.blastp.collinearity.kb.png")

### 基因组共线性展示
&ensp;&ensp;&ensp;&ensp;[circle_all]  
&ensp;&ensp;&ensp;&ensp;lens = lens file# 染色体文件  
&ensp;&ensp;&ensp;&ensp;gff = gff file# 基因注释文件  
&ensp;&ensp;&ensp;&ensp;species_list = Genome name# 物种列表例如：ath,vvi  
&ensp;&ensp;&ensp;&ensp;blockfile = block file# 共线性文件  
&ensp;&ensp;&ensp;&ensp;radius = 0.3# 半径  
&ensp;&ensp;&ensp;&ensp;savefile = save file (*.png, *.pdf)# 共线性展示图  
例：![test:拟南芥基因组共线性展示](https://images.gitee.com/uploads/images/2021/0717/135043_4acd1719_8074509.png "test.collinearity.circle.all.png")

### 局部基因共线性展示
&ensp;&ensp;&ensp;&ensp;[line]  
&ensp;&ensp;&ensp;&ensp;pairs_file = pairs file# 共线性文件  
&ensp;&ensp;&ensp;&ensp;gff1 =  gff1 file# 基因注释文件  
&ensp;&ensp;&ensp;&ensp;gff2 =  gff2 file# 基因注释文件  
&ensp;&ensp;&ensp;&ensp;lens1 = lens1 file# 染色体文件  
&ensp;&ensp;&ensp;&ensp;lens2 = lens2 file# 染色体文件  
&ensp;&ensp;&ensp;&ensp;chr1_name =  chr1 name# 展示染色体1  
&ensp;&ensp;&ensp;&ensp;chr2_name =  chr2 name# 展示染色体2  
&ensp;&ensp;&ensp;&ensp;savefile = savefile(.png,.pdf)# 局部共线性展示结果  
例：![test:拟南芥2号和3号染色体间共线性展示](https://images.gitee.com/uploads/images/2021/0717/131629_dd9e4e79_8074509.png "test.line.png")


## 家族展示
### 基因家族查找
&ensp;&ensp;&ensp;&ensp;[hmmer]  
&ensp;&ensp;&ensp;&ensp;pep = pep file# 蛋白质文件  
&ensp;&ensp;&ensp;&ensp;cds = cds file# 核酸文件(可选)  
&ensp;&ensp;&ensp;&ensp;hmmmoldpath = hmm file# 模型地址  
&ensp;&ensp;&ensp;&ensp;format_conversion = Fales# 格式转换  
&ensp;&ensp;&ensp;&ensp;comparison = muscle/clustal# 比对软件  
&ensp;&ensp;&ensp;&ensp;e_value1 = value1# 初筛阈值  
&ensp;&ensp;&ensp;&ensp;e_value2 = value2# 复筛阈值  

### 家族关系展示
&ensp;&ensp;&ensp;&ensp;[circle_family]  
&ensp;&ensp;&ensp;&ensp;lens1 = lens1 file  
&ensp;&ensp;&ensp;&ensp;lens2 = lens2 file  
&ensp;&ensp;&ensp;&ensp;gff1 = gff1 file  
&ensp;&ensp;&ensp;&ensp;gff2 = gff2 file  
&ensp;&ensp;&ensp;&ensp;species1 = species1 name  
&ensp;&ensp;&ensp;&ensp;species2 = species1 name  
&ensp;&ensp;&ensp;&ensp;genepairs = block file  
&ensp;&ensp;&ensp;&ensp;genepairsfile_type = MCScanX  
&ensp;&ensp;&ensp;&ensp;family_list = family file  
&ensp;&ensp;&ensp;&ensp;radius = 0.3  
&ensp;&ensp;&ensp;&ensp;block = 6# 当家族genepairfile为blast时block长度无意义  
&ensp;&ensp;&ensp;&ensp;savefile = save file (*.png, *.pdf)  
例：![拟南芥部分家族关系](https://images.gitee.com/uploads/images/2021/0918/114456_e4569355_8074509.png "ath.blast.cf.png")
![水稻NBS家族分布](https://images.gitee.com/uploads/images/2021/1011/165044_4086989e_8074509.png "osa_osa.wgdi.collinearity.wf.png")

### 结构域分布情况
&ensp;&ensp;&ensp;&ensp;[screen]  
&ensp;&ensp;&ensp;&ensp;domainpath = domain file path# 结构域搜索结果储存路径  
&ensp;&ensp;&ensp;&ensp;lens = lens file# 染色体文件  
&ensp;&ensp;&ensp;&ensp;gff = gff file# 基因注释文件  
&ensp;&ensp;&ensp;&ensp;chrolist = Genome name# 物种列表例如：ath,vvi  
&ensp;&ensp;&ensp;&ensp;series = 25# 串联数  
&ensp;&ensp;&ensp;&ensp;outpath = out file# 可视化结果储存路径  
例：![test:拟南芥1号染色体部分结构域分布情况](https://images.gitee.com/uploads/images/2021/0717/131138_a2661e2b_8074509.png "ath1.png")

### 文件格式化
&ensp;&ensp;&ensp;&ensp;[typing]  
&ensp;&ensp;&ensp;&ensp;domainpath = domain file# 结构域搜索结果储存路径  
&ensp;&ensp;&ensp;&ensp;domainlist = Genome name# 需要格式化的结构域搜索结果  
&ensp;&ensp;&ensp;&ensp;savefile = out file# 保存文件  

### 内卷型的tandem可视化
&ensp;&ensp;&ensp;&ensp;[inner]  
&ensp;&ensp;&ensp;&ensp;lens = lens file# 染色体文件  
&ensp;&ensp;&ensp;&ensp;gff = gff file# 基因注释文件  
&ensp;&ensp;&ensp;&ensp;chrolist = Genome name# 物种列表例如：ath,vvi  
&ensp;&ensp;&ensp;&ensp;ks = ks file# ks文件  
&ensp;&ensp;&ensp;&ensp;genefamily = famliy file# typing格式化的基因家族文件  
&ensp;&ensp;&ensp;&ensp;Ks_concern = 0,1.5# ks区间  
&ensp;&ensp;&ensp;&ensp;radius = 0.3# 半径  
&ensp;&ensp;&ensp;&ensp;space = 0.005# tandem基因间距  
&ensp;&ensp;&ensp;&ensp;clusters = True# ks区间内分段展示  
&ensp;&ensp;&ensp;&ensp;peripheral = False# 外围ks展示  
&ensp;&ensp;&ensp;&ensp;savecsv = outer file (*.csv)# 筛选结果储存文件  
&ensp;&ensp;&ensp;&ensp;savefile = save file (*.png, *.pdf)# 筛选结果可视化  
例：![test:拟南芥部分结构域tandem展示](https://images.gitee.com/uploads/images/2021/0717/132312_93086faa_8074509.png "test.inner.png")
![水稻内卷的家族分布图](https://images.gitee.com/uploads/images/2021/1011/165242_80621c7a_8074509.png "osa_osa.collinearity.ks.i.png")

### 放射型的tandem可视化
&ensp;&ensp;&ensp;&ensp;[outer]  
&ensp;&ensp;&ensp;&ensp;lens = lens file# 染色体文件  
&ensp;&ensp;&ensp;&ensp;gff = gff file# 基因注释文件  
&ensp;&ensp;&ensp;&ensp;chrolist = Genome name# 物种列表例如：ath,vvi  
&ensp;&ensp;&ensp;&ensp;ks = ks file# ks文件  
&ensp;&ensp;&ensp;&ensp;genefamily = famliy file# typing格式化的基因家族文件  
&ensp;&ensp;&ensp;&ensp;Ks_concern = 0,0.15# ks区间  
&ensp;&ensp;&ensp;&ensp;radius = 0.3# 半径  
&ensp;&ensp;&ensp;&ensp;space = 0.005# tandem基因间距  
&ensp;&ensp;&ensp;&ensp;clusters = True# ks区间内分段展示  
&ensp;&ensp;&ensp;&ensp;peripheral = False# 外围ks展示  
&ensp;&ensp;&ensp;&ensp;savecsv = outer file (*.csv)# 筛选结果储存文件  
&ensp;&ensp;&ensp;&ensp;savefile = save file (*.png, *.pdf)# 筛选结果可视化  
例：![test:拟南芥部分结构域tandem展示](https://images.gitee.com/uploads/images/2021/0717/135152_5bdfd054_8074509.png "test.outer.png")


### 局部tandem展示
&ensp;&ensp;&ensp;&ensp;[part]  
&ensp;&ensp;&ensp;&ensp;lens = lens file# 染色体文件  
&ensp;&ensp;&ensp;&ensp;gff = gff file# 基因注释文件  
&ensp;&ensp;&ensp;&ensp;chro_name = chro name# 展示染色体  
&ensp;&ensp;&ensp;&ensp;genefamily = genefamily file# typing格式化的基因家族文件  
&ensp;&ensp;&ensp;&ensp;pairs_file = pairs file# inner/outer筛选结果储存文件  
&ensp;&ensp;&ensp;&ensp;interval = 0,9000# 展示基因范围  
&ensp;&ensp;&ensp;&ensp;space = 0.005# tandem基因间距  
&ensp;&ensp;&ensp;&ensp;clusters = True# ks区间内分段展示  
&ensp;&ensp;&ensp;&ensp;savefile = save file (*.png, *.pdf)# 筛选结果可视化  
例：![test:拟南芥2号染色体部分结构域的tandem展示](https://images.gitee.com/uploads/images/2021/0717/133034_2af62f7f_8074509.png "test.part.png")
![水稻局部家族分布](https://images.gitee.com/uploads/images/2021/1011/165322_4ec10895_8074509.png "outer.out.csv.p.png")
