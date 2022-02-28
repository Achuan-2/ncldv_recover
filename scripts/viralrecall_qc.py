
#!/usr/bin/env python3

"""
viralrecall 质检标准
* 剔除污染的contig：只保留score > 0 & num_viralhits >3 & contig_length > 5000 过滤 bin.fa 文件中的contig
* 剔除污染contig后，只保留总大小超过 100 kb & mean score >1 & marker gene hit有3个及以上或要有mcp基因的bin
"""
import pandas as pd
import os.path
import argparse


def main(project, output):
    infer_NCLDV_num = 0
    table = []
    bin_dirs = traversalDir_FirstDir(f"{project}")

    for bin_dir in bin_dirs:
        #分箱水平
        bin_path = f"{project}/{bin_dir}/{bin_dir}"
        summary_table = f"{bin_path}.summary.tsv"
        origin_fasta = f"{bin_path}.fa"

        output_path = f"{output}/bin_dir"
        df = pd.read_table(summary_table)
        filter_df = df[(df["score"] > 0) & (df["num_viralhits"]
                                            >= 3) & (df["contig_length"] >= 5000)]
        filter_contig = len(df)-len(filter_df)
        filter_id = filter_df["replicon"]
        fasta_dict = read_fasta(origin_fasta)
        filtered_dict = {key: value for key,
                         value in fasta_dict.items() if key in list(filter_id)}
        # 计算总长度和平均分
        length_sum = 0
        for seq in filtered_dict.values():
            length_sum += len(seq)
        mean_score = filter_df["score"].mean()
        # 筛选bin：bin 的总长度要超过100kb且平均分要高于1才可以通过

        maker_set = set(filter_df["markerhits"])
        if '-' in maker_set:
            maker_set.remove('-')
        marker_list = []
        for i in maker_set:
            marker = i.split(':')[0]
            if marker not in marker_list:
                marker_list.append(marker)
        flag = ("mcp" in marker_list) | (len(marker_list) >= 3)
        if length_sum > 100000 and mean_score > 1 and flag:
            print("1")
            infer_NCLDV_num += 1
            row = {'bin': bin_dir, 'before_contig_num': df.shape[0], 'before_length': df["contig_length"].sum(), 'before_mean_score': df["score"].mean(),
                   'filter_contig_num': filter_contig, 'contig_num': filter_df.shape[0], 'length': length_sum, 'mean_score': mean_score}
            table.append(row)
            # 如果认为是NCLDV就生成得分pdf、统计tsv和过滤后的fasta序列
            mkdir(output_path)
            os.system(
                f"ln -s {bin_path}.pdf {output_path}/{bin_dir}.pdf ")
            filter_df.to_csv(
                f"{output_path}/{bin_dir}.summray_filter.tsv", sep="\t", index=False)
            write_fasta(f"{output_path}/{bin_dir}_filtered.fa", filtered_dict)

    print(f"共筛选出{infer_NCLDV_num}个NCLDV!")
    infer_NCLDV_df = pd.DataFrame(table)
    if not infer_NCLDV_df.empty:
        infer_NCLDV_df = infer_NCLDV_df.sort_values(
            by="mean_score", ascending=False)
        infer_NCLDV_df.to_csv(
            f"{output}/infer.tsv", sep="\t", index=False)


#定义一个函数，path为你的路径
def traversalDir_FirstDir(path):
    #定义一个列表，用来存储结果
    path_list = []
    #判断路径是否存在
    if (os.path.exists(path)):
        #获取该目录下的所有文件或文件夹目录
        files = os.listdir(path)
        for file in files:
            #得到该文件下所有目录的路径
            m = os.path.join(path, file)
            #判断该路径下是否是文件夹
            if (os.path.isdir(m)):
                h = os.path.split(m)
                path_list.append(h[1])
        return sorted(path_list)


def read_fasta(filename):
    dict = {}
    with open(filename, 'r') as fasta_f:
        for line in fasta_f:
            if line.startswith('>'):
                name = line.rstrip()[1:]
                dict[name] = ''
            else:
                dict[name] += line.rstrip()  # 读取整个fasta文件构成字典
    return dict


def write_fasta(output, fasta_dict):
    with open(output, 'w') as fasta_f:
        for key in fasta_dict.keys():  # 选取包含所需名称的基因名和序列
            fasta_f.write(">"+key + '\n')
            fasta_f.write(fasta_dict[key] + '\n')


def mkdir(path):
	folder = os.path.exists(path)
	if not folder:  # 判断是否存在文件夹如果不存在则创建为文件夹
		os.makedirs(path)  # makedirs 创建文件时如果路径不存在会创建这个路径




if __name__=="__main__":
    # 0. 读取参数
    # 创建一个ArgumentParser对象，以存储实参信息
    description = """
    viralrecall 质检标准
    * 剔除污染的contig：只保留score > 0 & num_viralhits >3 & contig_length > 5000 过滤 bin.fa 文件中的contig
    * 剔除污染contig后，只保留总大小超过 100 kb & mean score >1 & marker gene hit有4个及以上或要有mcp基因的bin
    """
    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawTextHelpFormatter)
    # 方法add_argument()添加要解析的命令内容
    parser.add_argument(
        '--input', '-i', type=str, help="input_dir: viralrecall result dir", required=True)
    parser.add_argument(
        '--output', '-o', type=str, help="output_dir: filtered output dir", required=True)
    args = parser.parse_args()  # 读入输入的参数，生成一个列表args
    work_dir = args.input  
    output = args.output
    # 1.主函数
    main(work_dir,output)
