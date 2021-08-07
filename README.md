# famCircle 基因家族分析，tandem可视化分析

## 分析流程
![分析流程图](https://images.gitee.com/uploads/images/2021/0717/124804_dc342bac_8074509.png "famCircle.png")
基因组下载以及处理
	路线：
		前期处理序列比对、共线性扫描、ks计算
		hmmer 隐马尔可夫模型搜索筛选
		screen tandem数量分部
			断开处：
				pfam结构域筛选# http://pfam.xfam.org/
				α/β筛选# http://cb.csail.mit.edu/cb/paircoil2/
		主程序选择：
			inner
				typing 数据格式化
				inner绘图
			outer
				typing 数据格式化
				outer绘图
		part局部tandem绘制
        其他可视化分析工具：
		ks展示
			全基因组 ks展示
			block ks展示
		circle全基因组共线性展示
		line局部共线性绘制


# 配置文件解释
## 基因家族查找
[hmmer]
pep = pep file# 蛋白质文件
cds = cds file# 核酸文件(可选)
hmmmoldpath = hmm file# 模型地址
format_conversion = Fales# 格式转换
comparison = muscle/clustal# 比对软件
e_value1 = value1# 初筛阈值
e_value2 = value2# 复筛阈值

## 结构域筛选
## α/β筛选

## 结构域分布情况
[screen]
domainpath = domain file path# 结构域搜索结果储存路径
lens = lens file# 染色体文件
gff = gff file# 基因注释文件
chrolist = Genome name# 物种列表例如：ath,vvi
series = 25# 串联数
outpath = out file# 可视化结果储存路径
例：![test:拟南芥1号染色体部分结构域分布情况](https://images.gitee.com/uploads/images/2021/0717/131138_a2661e2b_8074509.png "ath1.png")

## 文件格式化
[typing]
domainpath = domain file# 结构域搜索结果储存路径
domainlist = Genome name# 需要格式化的结构域搜索结果
savefile = out file# 保存文件

## 内卷型的tandem可视化
[inner]
lens = lens file# 染色体文件
gff = gff file# 基因注释文件
chrolist = Genome name# 物种列表例如：ath,vvi
ks = ks file# ks文件
genefamily = famliy file# typing格式化的基因家族文件
Ks_concern = 0,1.5# ks区间
radius = 0.3# 半径
space = 0.005# tandem基因间距
clusters = True# ks区间内分段展示
peripheral = False# 外围ks展示
savecsv = outer file (*.csv)# 筛选结果储存文件
savefile = save file (*.png, *.pdf)# 筛选结果可视化
例：![test:拟南芥部分结构域tandem展示](https://images.gitee.com/uploads/images/2021/0717/132312_93086faa_8074509.png "test.inner.png")

## 放射型的tandem可视化
[outer]
lens = lens file# 染色体文件
gff = gff file# 基因注释文件
chrolist = Genome name# 物种列表例如：ath,vvi
ks = ks file# ks文件
genefamily = famliy file# typing格式化的基因家族文件
Ks_concern = 0,0.15# ks区间
radius = 0.3# 半径
space = 0.005# tandem基因间距
clusters = True# ks区间内分段展示
peripheral = False# 外围ks展示
savecsv = outer file (*.csv)# 筛选结果储存文件
savefile = save file (*.png, *.pdf)# 筛选结果可视化
例：![test:拟南芥部分结构域tandem展示](https://images.gitee.com/uploads/images/2021/0717/135152_5bdfd054_8074509.png "test.outer.png")

## 局部tandem展示
[part]
lens = lens file# 染色体文件
gff = gff file# 基因注释文件
chro_name = chro name# 展示染色体
genefamily = genefamily file# typing格式化的基因家族文件
pairs_file = pairs file# inner/outer筛选结果储存文件
interval = 0,9000# 展示基因范围
space = 0.005# tandem基因间距
clusters = True# ks区间内分段展示
savefile = save file (*.png, *.pdf)# 筛选结果可视化
例：![test:拟南芥2号染色体部分结构域的tandem展示](https://images.gitee.com/uploads/images/2021/0717/133034_2af62f7f_8074509.png "test.part.png")



## 计算ks
[Ks]
cds_file = 	cds file# 基因组cds文件
pep_file = 	pep file# 基因组蛋白质文件(可选)
align_software = muscle# 比对软件
pairs_file = gene pairs file# 共线性文件
ks_file = ks result# 输出文件

## 全基因组ks数据分布
[Ks_allocation]
ks = ks file# ks文件
species_list = specise name# 物种列表例如：ath,vvi
area = 0,2.5# ks绘制区间
model = YN00/NG86# 绘制ks来源
savefile = save file (*.png, *.pdf)# 输出图片
例：![test:拟南芥ks数据分布](https://images.gitee.com/uploads/images/2021/0717/125414_357dd5c8_8074509.png "test.collinearity.ks.all.png")

## 基因块ks数据分布图
[Ks_block]
species_list = specise name# 物种列表例如：ath,vvi
ks = ks file# ks文件
area = 0,2.5# ks绘制区间
model = YN00/NG86# 绘制ks来源
blockfile = block file# 共线性文件
blocklength = 6# block长度控制参数
pvalue = 1e-5# block阈值控制参数
savecsv = save csv# 筛选结果文件
savefile = save file (*.png, *.pdf)# 筛选结果可视化
例：![test:拟南芥基因块ks数据分布block长度为10+](https://images.gitee.com/uploads/images/2021/0717/125847_77aaf58a_8074509.png "test.collinearity.ks.block.png")

## 基因组共线性展示
[circle_all]
lens = lens file# 染色体文件
gff = gff file# 基因注释文件
species_list = Genome name# 物种列表例如：ath,vvi
blockfile = block file# 共线性文件
radius = 0.3# 半径
savefile = save file (*.png, *.pdf)# 共线性展示图
例：![test:拟南芥基因组共线性展示](https://images.gitee.com/uploads/images/2021/0717/135043_4acd1719_8074509.png "test.collinearity.circle.all.png")

## 局部基因共线性展示
[line]
pairs_file = pairs file# 共线性文件
gff1 =  gff1 file# 基因注释文件
gff2 =  gff2 file# 基因注释文件
lens1 = lens1 file# 染色体文件
lens2 = lens2 file# 染色体文件
chr1_name =  chr1 name# 展示染色体1
chr2_name =  chr2 name# 展示染色体2
savefile = savefile(.png,.pdf)# 局部共线性展示结果
例：![test:拟南芥2号和3号染色体间共线性展示](https://images.gitee.com/uploads/images/2021/0717/131629_dd9e4e79_8074509.png "test.line.png")

## 结合ks的共线性展示
[circle]
lens = lens file# 染色体文件
gff = gff file# 基因注释文件
species_list = Genome name# 物种列表例如：ath,vvi
ks = ks file# ks文件
genepairs = genepairs file# 共线性文件
block = 6# block长度控制参数
radius = 0.45# 半径
savefile = save file (*.png, *.pdf)# ks共线性展示图
例：![test:拟南芥ks共线性可视化](https://images.gitee.com/uploads/images/2021/0717/135435_0c4f7f81_8074509.png "test.collinearity.circle.png")

