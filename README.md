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
ks = ks file
ks_concern = min,max
step_size = 0.05
model = YN00
savefile = save file (*.png, *.pdf)

## 基因家族查找
run hmmer 
pep = pep file# 蛋白质文件
hmmmoldpath = hmm file# 模型地址
format_conversion = Fales# 格式转换
comparison = clustal# 比对软件
e_value1 = value1
e_value2 = value2
structure_field = False# 结构域

## 内卷型的tandem可视化
run inner
lens = lens file
gff = gff file
chrolist = Genome name
ks = ks file
genefamily = famliy file
Ks_concern = 0.1,0.2,1.4,1.5# ks分割参数可调
peripheral = True# ks显示范围参数
bridge = 0.05# 最远连接
savecsv = outer file (*.csv)
savefile = save file (*.png, *.pdf)

## 共线性可视化
run circle
lens = lens file
gff = gff file
chrolist = Genome name
ks = ks file
genepairs = genepairs file
Ks_concern = 0,0.15
peripheral = False
bridge = 1
savefile = save file (*.png, *.pdf)

## 放射型的tandem可视化
run outer
lens = lens file
gff = gff file
chrolist = Genome name
ks = ks file
genepairs = genepairs file
Ks_concern = 0,0.15# ks分割参数可调
clusters = 100000# 集团参数可调
peripheral = False# ks范围显示
bridge = 1# 跨区域连线参数可调
savecsv = outer file (*.csv)
savefile = save file (*.png, *.pdf)
