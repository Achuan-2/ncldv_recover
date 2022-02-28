#!/usr/bin/env python3

"""help
input：work_dir绝对路径
脚本： /mnt/raid7/Dachuang/Achuan/scripts/viralrecall_filter.py
功能：
viralrecall使用-db ""marker"后筛选hmm，

根据平均分cutoff值（设置为1）和markhit数目超过或等于3或包含mcp基因的bin，才会筛选进infer_NCLDV.tsv中;如果没有符合条件的bin，则没有infer_NCLDV.tsv这个文件
usage:
# 一个文件夹
python /mnt/raid7/Dachuang/Achuan/scripts/viralrecall_filter.py --i /mnt/raid7/Dachuang/Achuan/02_is_NCLDV/viralrecall/bins_out 

## 多个文件夹
PROJECT=/mnt/raid7/Dachuang/Achuan/02_is_NCLDV/viralrecall_result
DATA=data15
ls ${PROJECT}/${DATA} |while read id;do python viralrecall_summary.py --i ${PROJECT}/${DATA}/${ID} ;done
"""

# 0. 读取参数
# 创建一个ArgumentParser对象，以存储实参信息

import os
import pandas as pd
import argparse


def main():
    # 0. 读取参数
    # 创建一个ArgumentParser对象，以存储实参信息
    description = """
    功能：
    viralrecall使用-db ""marker"后筛选hmm，

    根据平均分cutoff值（设置为1）和markhit数目超过3或包含mcp基因的bin，才会筛选进infer_NCLDV.tsv中;如果没有符合条件的bin，则没有infer_NCLDV.tsv这个文件
    """
    parser = argparse.ArgumentParser(
        description=description, formatter_class=argparse.RawTextHelpFormatter)
    # 方法add_argument()添加要解析的命令内容
    parser.add_argument(
        '--input', '-i', type=str, help="input_dir: viralrecall result dir", required=True)
    args = parser.parse_args()  # 读入输入的参数，生成一个列表args
    work_dir = args.input  # 接着对参数的任何操作，调用命名为xxx的参数方式为args.xxx

    # 1. 根据markhit数目超过3或包含mcp基因进行bin筛选
    filter(work_dir)


def filter(work_dir):
    summary_df=pd.read_table(work_dir+"/summary.tsv", sep="\t")
    summary_df["markerhits"] = summary_df["markerhits"]
    flag = (summary_df["markerhits"].str.contains("mcp")) | (summary_df["markerhits_num"] >=3)


    if sum(flag):
        infer_NCDLV = summary_df[flag]
        infer_NCDLV.to_csv(work_dir+"/marker_filtered_bin.tsv", sep="\t", index=False)

if __name__ == '__main__':
    main()
