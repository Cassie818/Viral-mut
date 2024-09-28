#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
clean_fasta.py

这个脚本用于处理FASTA格式的文件，去除所有空格和换行符，并将所有核苷酸字符转换为大写。
"""

import argparse
import sys


def clean_fasta(input_file, output_file):
    """
    处理FASTA文件，去除空格和换行符，并将核苷酸字符转换为大写。

    参数:
        input_file (str): 输入的FASTA文件路径。
        output_file (str): 输出的处理后FASTA文件路径。
    """
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
            sequence = []
            header = None
            line_number = 0  # 用于调试
            for line in infile:
                line_number += 1
                line = line.strip()
                if not line:
                    continue  # 跳过空行
                if line.startswith('>'):
                    if header:
                        # 写入前一个序列
                        outfile.write(header + '\n')
                        outfile.write(''.join(sequence) + '\n')
                        print(f"已处理序列: {header}，长度: {len(''.join(sequence))}")
                        sequence = []
                    header = line
                else:
                    # 去除所有空格并转换为大写
                    cleaned_seq = line.replace(" ", "").upper()
                    sequence.append(cleaned_seq)
            # 写入最后一个序列
            if header and sequence:
                outfile.write(header + '\n')
                outfile.write(''.join(sequence) + '\n')
                print(f"已处理序列: {header}，长度: {len(''.join(sequence))}")
        print(f"\n处理完成！已生成文件：{output_file}")
    except FileNotFoundError:
        print(f"错误：无法找到文件 '{input_file}'。请检查文件路径。", file=sys.stderr)
    except Exception as e:
        print(f"发生错误：{e}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="处理FASTA文件，去除空格和换行符，并将核苷酸字符转换为大写。")
    parser.add_argument('-i', '--input', required=True, help="输入的FASTA文件路径")
    parser.add_argument('-o', '--output', required=True, help="输出的处理后FASTA文件路径")
    args = parser.parse_args()
    clean_fasta(args.input, args.output)


if __name__ == "__main__":
    main()
