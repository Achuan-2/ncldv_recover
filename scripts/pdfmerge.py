#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
本脚本用来合并viralrecall预测得分的pdf文件，输出的pdf文件按输入的pdf文件名生成书签
对于单个样本，使用示例如下：
python ${SCRIPTS}/pdfmerge.py -i "/mnt/raid7/Dachuang/Achuan/02_is_NCLDV/viralrecall/bins_out" -o "merged.pdf" -b False

示例说明：
要合并的pdf文件所在的路径： /mnt/raid7/Dachuang/Achuan/02_is_NCLDV/viralrecall/bins_out
合并后的pdf文件的输出文件名：merged.pdf
是否从pdf文件中导入书签的值：False

'''
import os
import sys
import codecs
from argparse import ArgumentParser, RawTextHelpFormatter
from PyPDF2 import PdfFileReader, PdfFileWriter, PdfFileMerger

def main():
    description = """
    本脚本用来合并viralrecall预测得分的pdf文件，输出的pdf文件按输入的pdf文件名生成书签
    对于单个样本，使用示例如下：
    python /mnt/raid7/Dachuang/Achuan/scripts/pdfmerge.py -i "/mnt/raid7/Dachuang/Achuan/02_is_NCLDV/viralrecall/bins_out" -o "merged.pdf" -b False

    示例说明：
    要合并的pdf文件所在的路径： /mnt/raid7/Dachuang/Achuan/02_is_NCLDV/viralrecall/bins_out
    合并后的pdf文件的输出文件名：merged.pdf
    是否从pdf文件中导入书签的值：False

    对于数据集下的所有样本
    PROJECT=/mnt/raid7/Dachuang/Achuan/02_is_NCLDV/viralrecall_result
    DATA=data15
    ls ${PROJECT}/${DATA} |while read id;do python3 /mnt/raid7/Dachuang/Achuan/scripts/pdfmerge.py --i ${PROJECT}/${DATA}/${ID}  -o "merged.pdf";done
    """

    # 添加程序帮助，程序帮助支持换行符号
    parser = ArgumentParser(description=description,
                            formatter_class=RawTextHelpFormatter)

    # 添加命令行选项

    parser.add_argument("-i", "--input",
                        dest="input_path",
                        default=".",
                        help="PDF文件所在目录")
    parser.add_argument("-o", "--output",
                        dest="output_filename",
                        default="merged.pdf",
                        help="合并PDF的输出文件名",
                        metavar="FILE")
    parser.add_argument("-b", "--bookmark",
                        dest="import_bookmarks",
                        default="False",
                        help="是否从pdf文件中导入书签，值可以是'True'或者'False'")

    args = parser.parse_args()
    try:
        mergefiles(args.input_path, args.output_filename,
                args.import_bookmarks)
    except:
        print('Error to merge pdf file:')
        print(sys.exc_info()[0],sys.exc_info()[1])

def getfilenames(filepath='', filelist_out=[], file_ext='.pdf'):
    # 遍历filepath下的所有文件，包括子目录下的文件, 获取pdf
    for fpath, dirs, fs in os.walk(filepath):
        for f in fs:
            fi_d = os.path.join(fpath, f)
            if os.path.splitext(fi_d)[1] == file_ext:
                filelist_out.append(fi_d)
            else:
                pass
    return filelist_out

def mergefiles(path, output_filename, import_bookmarks=False):
    # 遍历目录下的所有pdf将其合并输出到一个pdf文件中，输出的pdf文件默认带书签，书签名为之前的文件名
    # 默认情况下原始文件的书签不会导入，使用import_bookmarks=True可以将原文件所带的书签也导入到输出的pdf文件中
    merger = PdfFileMerger()
    filelist = getfilenames(filepath=path, file_ext='.pdf')
    out_filename = os.path.join(os.path.abspath(path), output_filename)
    if len(filelist) == 0:
        print("当前目录及子目录下不存在pdf文件")
        sys.exit()
    for filename in filelist:
        if filename == out_filename:
            print("已存在合并的pdf，将覆盖原来文件")
            continue
        f = codecs.open(filename, 'rb')
        file_rd = PdfFileReader(f)
        short_filename = os.path.basename(os.path.splitext(filename)[0])
        if file_rd.isEncrypted == True:
            print('不支持的加密文件：%s' % (filename))
            continue
        merger.append(file_rd, bookmark=short_filename,
                        import_bookmarks=import_bookmarks)
        print('合并文件：%s' % (filename))
        f.close()
    merger.write(out_filename)
    print('合并后的输出文件：%s' % (out_filename))
    merger.close()


if __name__ == "__main__":
    main()
