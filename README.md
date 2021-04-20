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

## 基因家族查找
run hmmer 
pep = pep file
hmmmold = hmm file
comparison = clustal# 比对软件
e_value1 = value1# e值筛选
e_value2 = value2
structure_field = True/False# 结构域筛选
newmold = hmm file# 本地生成隐马尔可夫模型
savealn = aln file (*.aln)# 生成的树文件
savefile = family pep
hmmlist = genefamily list (*.out)# 基因家族列表文件

## ks分布可视化
ks = ks file
ks_concern = min,max
step_size = 0.05
model = YN00
savefile = save file (*.png, *.pdf)

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
