# famCircle 基因家族分析，tandem可视化分析

## 分析流程

基因组下载以及处理
	路线1：
		ks展示
		hmmer搜索筛选
			断开处：
				pfam结构域筛选
				α/β筛选
		主程序选择：
			inner
			outer
	路线2：
		blast、共线性文件筛选
		主程序选择：
			inner
			outer

## ks分布可视化
run ks
ks = ks file
vertical = False
bins = 100
model = YN00/NG86
savefile = save file (*.png, *.pdf)

## 结构域分布情况
run screen
domainpath = domain file
lens = lens file
gff = gff file
chrolist = Genome name
series = 25
outpath = out file

## 基因家族查找
run hmmer 
pep = pep file# 蛋白质文件
hmmmoldpath = hmm file# 模型地址
format_conversion = Fales# 格式转换
comparison = muscle/clustal# 比对软件
e_value1 = value1
e_value2 = value2

## 结构域筛选

## α/β筛选

## 内卷型的tandem可视化
run inner
lens = lens file
gff = gff file
chrolist = Genome name
ks = ks file
genefamily = famliy file
Ks_concern = 0.1,0.2,1.4,1.5# ks分割参数可调
bridge = 0.05# 最远连接
radius = 0.3
peripheral = False
savecsv = outer file (*.csv)
savefile = save file (*.png, *.pdf)

## 共线性可视化
run circle
lens = lens file
gff = gff file
chrolist = Genome name
ks = ks file
genepairs = genepairs file
peripheral = False
block = 6
Ks_concern = 0,0.15
bridge = 1
radius = 0.3
savefile = save file (*.png, *.pdf)

## 放射型的tandem可视化
run outer
lens = lens file
gff = gff file
chrolist = Genome name
ks = ks file
genepairs = genepairs file
Ks_concern = 0,0.15# ks分割参数可调
series = 25# 串联数
clusters = None# 集团参数可调
peripheral = False
bridge = 1# 跨区域连线参数可调
radius = 0.3
savecsv = outer file (*.csv)
savefile = save file (*.png, *.pdf)
