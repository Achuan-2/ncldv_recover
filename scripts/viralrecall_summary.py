#!/usr/bin/env python3

"""help
input：work_dir绝对路径
脚本： /mnt/raid7/Dachuang/Achuan/scripts/viralrecall_filter.py
功能：
viralrecall基于MAG水平筛选NCLDV，
1. 先统计每个bin的平均得分，contig数目，总长度，最大contig长度，最小contig长度,makerhit数目
2. 根据平均分cutoff值（设置为1）和markhit数目超过3或包含mcp基因的bin，才会筛选进infer_NCLDV.tsv中;如果没有符合条件的bin，则没有infer_NCLDV.tsv这个文件

usage:
# 一个文件夹
python /mnt/raid7/Dachuang/Achuan/scripts/viralrecall_filter.py --i /mnt/raid7/Dachuang/Achuan/02_is_NCLDV/viralrecall/bins_out 

## 多个文件夹
PROJECT=/mnt/raid7/Dachuang/Achuan/02_is_NCLDV/viralrecall_result
DATA=data15
ls ${PROJECT}/${DATA} |while read id;do python viralrecall_summary.py --i ${PROJECT}/${DATA}/${ID} ;done
"""


import os
import pandas as pd
import argparse


def main():

    # 0. 读取参数
    # 创建一个ArgumentParser对象，以存储实参信息
    description = """
    功能：
    对viralrecall单个样本不同分箱的结果进行汇总
    统计每个bin的平均得分，contig数目，总长度，最大contig长度，最小contig长度,makerhit数目
    """
    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawTextHelpFormatter)
    # 方法add_argument()添加要解析的命令内容
    parser.add_argument(
        '--input', '-i', type=str, help="input_dir: viralrecall result dir", required=True)
    args = parser.parse_args()  # 读入输入的参数，生成一个列表args
    work_dir = args.input  # 接着对参数的任何操作，调用命名为xxx的参数方式为args.xxx
    # 1. 先统计每个bin的平均得分
    summary(work_dir)

def summary(work_dir):
    table_df = []
    for parent, dirnames, filenames in os.walk(work_dir):
        for dir in dirnames:
            table = pd.read_table(f"{work_dir}/{dir}/{dir}.summary.tsv")
            mean_score = table["score"].mean()
            sum_length = table["contig_length"].sum()
            contig_num = table.shape[0]
            max_contig = table["contig_length"].max()
            min_contig = table["contig_length"].min()
            maker_set = set(table["markerhits"])
            if '-' in maker_set:
                maker_set.remove('-')
            marker_list = []
            for i in maker_set:
                marker = i.split(':')[0]
                if marker not in marker_list:
                    marker_list.append(marker)
            row={"bin": dir, "mean_score": mean_score, "sum_length": sum_length,\
                        "contig_num": contig_num, "max_contig": max_contig, "min_contig": min_contig,"markerhits": marker_list,"markerhits_num":len(marker_list)}
            # print(row)
            table_df.append(row)
    summary_score = pd.DataFrame(table_df)
    summary_score.to_csv(work_dir+"/summary.tsv", sep="\t", index=False)



if __name__ =="__main__":
    main()
