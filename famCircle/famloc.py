import pandas as pd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import os
from famCircle.bez import *

class famloc():
    def __init__(self, options):
        # 输入文件路径
        self.lens_file = "she.lens"
        self.gff_file = "she.gff"
        self.gene_list_file = "gene_list.txt"
        self.dpi = 600
        self.font_size=10
        self.text_size = 1.3
        self.savefile="chromosome_map_dynamic_spacing.png"
        for k, v in options:
            setattr(self, str(k), v)
            print(k, ' = ', v)
        self.dpi = int(self.dpi)
        self.font_size = int(self.font_size)
        self.text_size = float(self.text_size)

    # 读取 lens 文件
    def read_lens_file(self,lens_file):
        lens_df = pd.read_csv(lens_file, sep="\t", header=None, names=["chromosome", "length", "gene_count"])
        return lens_df

    # 读取 gff 文件
    def read_gff_file(self,gff_file):
        gff_df = pd.read_csv(gff_file, sep="\t", header=None, 
                             names=["chromosome", "gene_name", "start", "end","p", "relative_position","old"])
        return gff_df

    # 读取基因 list 文件
    def read_gene_list(self,gene_list_file):
        gene_list_df = pd.read_csv(gene_list_file, sep="\t", header=None, names=["gene_name"])
        return set(gene_list_df["gene_name"])

    # 计算文本尺寸
    def calculate_text_size(self,text, font_size, dpi, font_prop):
        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
        fig, ax = plt.subplots(figsize=(10, 10), dpi=dpi)
        # 锁定比例为 1:1
        ax.set_aspect('equal')
        ax.set_ylim(0, 1.2)
        ax.set_ylim(0, 1.2)
        canvas = FigureCanvas(fig)
        text_obj = fig.text(0, 0, text, fontsize=font_size, fontproperties=font_prop)
        bbox = text_obj.get_window_extent(renderer=canvas.get_renderer())
        # plt.savefig('text.png', dpi=dpi)
        plt.close(fig)
        return bbox.width / dpi, bbox.height / dpi  # 返回宽和高（以英寸为单位）

    # 绘制染色体和标记基因
    def plot_chromosomes_with_genes(self,lens_df, gff_df, gene_list, font_size, dpi, output_file):
        path = get_path()
        font_path = os.path.join(path, 'example/arial.ttf')
        print(font_path)
        from matplotlib.font_manager import FontProperties
        font_prop = FontProperties(fname=font_path)

        fig, ax = plt.subplots(figsize=(10, 10), dpi=dpi)
        # 锁定比例为 1:1
        ax.set_aspect('equal')
        maxlength = lens_df["length"].max()

        # 设置图形范围
        ax.set_ylim(-0.1, 1.2)
        ax.set_xlim(0, 1.4)
        # 隐藏 X 轴和右边框
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xticklabels([])
        ax.spines['bottom'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        # 绘制左侧箭头
        arrow_x = 0.1  # 箭头位置
        ax.annotate(
            '', xy=(arrow_x, 0), xytext=(arrow_x, 1),
            arrowprops=dict(arrowstyle="<|-", color="black", lw=4)
        )
        # 绘制刻度
        num_ticks = 6  # 刻度数量
        tick_positions = [i * maxlength * 1.2 / (num_ticks - 1) for i in range(num_ticks)]
        tick_positionsorder = [i * 1.2 / (num_ticks - 1) for i in range(num_ticks)]
        # 使用科学计数法格式化刻度标签
        tick_labels = [f"{pos:.1e}" for pos in tick_positions]  # 科学计数法格式

        for tick, label in zip(tick_positionsorder[:-1], tick_labels[:-1]):
            # 短横线位置和宽度
            ax.hlines(y=tick, xmin=arrow_x, xmax=arrow_x+0.01, color="black", linewidth=2)
            ax.text(arrow_x - 0.1, tick, label, ha="left", va="center", fontsize=10, color="black", fontproperties=font_prop)
        ax.text(arrow_x - 0.1, 1, 'bp', ha="left", va="center", fontsize=10, color="black", fontproperties=font_prop)
        # 计算最长基因名字所占的宽度和高度
        text_width, text_height = self.calculate_text_size("A" * max(len(gene) for gene in gene_list), font_size, dpi,font_prop)
        # 图形的实际尺寸
        fig_width, fig_height = fig.get_size_inches()
        # print(fig_width, fig_height,text_width, text_height)
        x_range = 1.2  # x 轴范围（染色体数量 * 每个染色体的宽度）
        y_range = 1.2  # y 轴范围（染色体最大长度的 1.2 倍）
        text_width_in_coords = (text_width * x_range) / fig_width
        text_height_in_coords = (text_height * y_range) / fig_height
        # 计算染色体之间的间距
        spacing = text_width_in_coords * self.text_size  # 比最长基因名字宽度多 50%
        x_pos = [0.2+i * spacing for i in range(len(lens_df))]
        # print(spacing,x_pos)

        label_positions = {}  # 存储基因名字的最终位置以避免重叠
        for i, row in lens_df.iterrows():
            x = x_pos[i]
            chrom = row["chromosome"]
            length = row["length"]/maxlength
            # print(chrom,length,x)
            ax.vlines(x=x, ymin=0, ymax=length, linewidth=5, color="gray")
            ax.text(x, -0.05, chrom, ha="center", va="center", fontsize=10, color="black", fontproperties=font_prop)
            chrom_genes = gff_df[gff_df["chromosome"] == chrom]
            gene_positions = []
            for _, gene_row in chrom_genes.iterrows():
                if gene_row["gene_name"] in gene_list:
                    gene_start = gene_row["start"]
                    gene_name = gene_row["gene_name"]
                    gene_positions.append((gene_start, gene_name))
            # 对基因位置排序以便于调整标签位置
            gene_positions.sort()
            last_label_y = -float('inf')  # 记录上一个标签的 Y 坐标
            for gene_start, gene_name in gene_positions:
                gene_start = gene_start/maxlength
                # 确保基因标签之间的距离足够大
                label_y = max(gene_start, last_label_y + (text_height_in_coords * 1.1))
                last_label_y = label_y
                label_positions[gene_name] = label_y
                # 绘制基因位置为点
                ax.plot(x, gene_start, 'ro')  # 红点表示基因位置
                # 绘制从基因位置到标签的线段
                ax.plot([x, x + text_width_in_coords * 0.1], [gene_start, label_y], 'r-', linewidth=0.8)
                # 绘制标签
                ax.text(x + text_width_in_coords * 0.1, label_y, gene_name, fontsize=font_size, color="red", va="center")
        plt.tight_layout()
        plt.savefig(output_file, dpi=dpi, bbox_inches='tight')

    # 主程序
    def run(self):
        # 读取文件
        lens_df = self.read_lens_file(self.lens_file)
        # print(lens_df)
        gff_df = self.read_gff_file(self.gff_file)
        # print(gff_df)
        gene_list = self.read_gene_list(self.gene_list_file)

        # 可视化基因位置
        self.plot_chromosomes_with_genes(lens_df, gff_df, gene_list,self.font_size,self.dpi,self.savefile)
