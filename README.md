# famCircle 基因家族分析，tandem可视化分析

## 基因家族查找
run hmmer 
pep = pep file
hmmmold = hmm file
newmold = hmm file
comparison = clustal
savealn = aln file (*.aln)
hmmlist = genefamily list (*.out)
savefile = family pep

## ks分布可视化
run lookKs
ks = ks file
Ks_concern = 0,0.15
peripheral = False
bridge = 1
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
