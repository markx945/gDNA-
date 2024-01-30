#!/bin/bash

# 显示帮助信息
function show_help {
    echo "Usage: $0 -i input_path [-o output_file]"
    echo
    echo "-i input_path    The path where to look for directories."
    echo "-o output_file   Optional. The file where to save the directory names."
    echo "                 If not provided, output will be printed to the console."
    echo
    echo "Example: $0 -i /path/to/directory -o /path/to/output.txt"
    echo
    exit 1
}

# 初始化输入和输出变量
input_path=""
output_file=""

# 解析命令行选项和参数
while getopts "hi:o:" opt; do
    case "$opt" in
    h)
        show_help
        ;;
    i)
        input_path=$OPTARG
        ;;
    o)
        output_file=$OPTARG
        ;;
    *)
        show_help
        ;;
    esac
done

# 检查是否提供了输入路径
if [ -z "$input_path" ]; then
    echo "Error: No input path provided."
    show_help
fi

# 读取路径下的所有文件夹并保存到数组中
directories=()
while IFS=  read -r -d $'\0'; do
    dir_name=$(basename "$REPLY")
    directories+=("$dir_name")
done < <(find "$input_path" -mindepth 1 -maxdepth 1 -type d -print0)

# # 输出文件夹列表
# if [ -z "$output_file" ]; then
#     # 没有提供输出文件，打印到控制台
#     echo "Directories in '$input_path':"
#     printf "%s\n" "${directories[@]}"
# else
#     # 提供了输出文件，写入文件
#     printf "%s\n" "${directories[@]}" > "$output_file"
#     echo "Directories have been saved to '$output_file'"
# fi



