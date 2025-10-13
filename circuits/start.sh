#!/bin/bash

# 检查参数是否正确
if [ $# -ne 2 ]; then
    echo "使用方法: $0 <程序路径> <目标文件夹>"
    echo "示例: $0 ./main ./test_files"
    exit 1
fi

PROGRAM="$1"
TARGET_DIR="$2"

# 检查程序是否存在且可执行
if [ ! -x "$PROGRAM" ]; then
    echo "错误: 程序 $PROGRAM 不存在或不可执行"
    exit 1
fi

# 检查目标文件夹是否存在
if [ ! -d "$TARGET_DIR" ]; then
    echo "错误: 文件夹 $TARGET_DIR 不存在"
    exit 1
fi

# 创建结果日志文件（所有输出将保存到这里）
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="benchmark_results_$TIMESTAMP.log"
echo "基准测试开始于: $(date)" > "$LOG_FILE"
echo "程序路径: $PROGRAM" >> "$LOG_FILE"
echo "目标文件夹: $TARGET_DIR" >> "$LOG_FILE"
echo "========================================" >> "$LOG_FILE"
echo "" >> "$LOG_FILE"

# 统计文件数量
FILE_COUNT=$(find "$TARGET_DIR" -type f | wc -l)
echo "发现 $FILE_COUNT 个文件，开始基准测试..."
echo "所有输出将保存到: $LOG_FILE"

# 遍历文件夹中的所有文件
FILE_NUM=1
find "$TARGET_DIR" -type f | while read -r file; do
    echo "正在处理文件 ($FILE_NUM/$FILE_COUNT): $file"

    # 记录开始时间
    START_TIME=$(date +%s%N)

    # 记录当前处理的文件信息到日志
    echo "----------------------------------------" >> "$LOG_FILE"
    echo "文件: $file" >> "$LOG_FILE"
    echo "运行时间: $(date)" >> "$LOG_FILE"
    echo "程序输出:" >> "$LOG_FILE"

    # 运行程序，将所有输出（stdout和stderr）都重定向到日志文件
    "$PROGRAM" "$file" >> "$LOG_FILE" 2>&1

    # 计算运行时间
    END_TIME=$(date +%s%N)
    DURATION=$(( (END_TIME - START_TIME) / 1000000 ))  # 转换为毫秒

    echo "运行时间: $DURATION 毫秒" >> "$LOG_FILE"
    echo "----------------------------------------" >> "$LOG_FILE"
    echo "" >> "$LOG_FILE"

    echo "完成处理，耗时: $DURATION 毫秒"
    FILE_NUM=$((FILE_NUM + 1))
done

echo "========================================" >> "$LOG_FILE"
echo "基准测试完成于: $(date)" >> "$LOG_FILE"
echo "所有结果已保存到: $LOG_FILE"

exit 0
