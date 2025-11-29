#!/bin/bash

# 检查参数是否正确
if [ $# -ne 2 ]; then
    echo "使用方法: $0 <程序路径> <目标文件夹>"
    echo "示例: $0 ./main ./test_files"
    exit 1
fi

PROGRAM="$1"
TARGET_DIR="$2"
MAX_MEM_MB=$((32 * 1024))  # 32GB = 32768MB

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

# 检查必要工具是否可用
if ! command -v /usr/bin/time &> /dev/null; then
    echo "错误: 未找到 /usr/bin/time 工具，无法监控内存占用"
    exit 1
fi

if ! command -v bc &> /dev/null; then
    echo "错误: 未找到 bc 工具，无法进行内存计算"
    exit 1
fi

# 创建结果日志文件
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
LOG_FILE="benchmark_results_$TIMESTAMP.log"
echo "基准测试开始于: $(date)" > "$LOG_FILE"
echo "程序路径: $PROGRAM" >> "$LOG_FILE"
echo "目标文件夹: $TARGET_DIR" >> "$LOG_FILE"
echo "内存限制: ${MAX_MEM_MB}MB (32GB)" >> "$LOG_FILE"
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

    # 记录当前处理的文件信息到日志
    echo "----------------------------------------" >> "$LOG_FILE"
    echo "文件: $file" >> "$LOG_FILE"
    echo "运行时间: $(date)" >> "$LOG_FILE"
    echo "程序输出:" >> "$LOG_FILE"

    # 创建临时文件存储程序输出
    OUTPUT_TMP=$(mktemp)
    TIME_TMP=$(mktemp)

    # 记录开始时间
    START_TIME=$(date +%s%N)

    # 启动程序并监控内存，使用子shell进行内存限制检查
    (
        # 监控子进程内存占用
        while true; do
            # 获取子进程ID和内存占用（KB）
            PID=$!
            if [ -z "$PID" ] || ! ps -p "$PID" >/dev/null; then
                break
            fi

            # 获取当前内存占用（KB）并转换为MB
            MEM_KB=$(ps -p "$PID" -o rss=)
            MEM_MB=$(echo "scale=2; $MEM_KB / 1024" | bc)

            # 检查是否超过内存限制
            if (( $(echo "$MEM_MB > $MAX_MEM_MB" | bc -l) )); then
                kill -9 "$PID" >/dev/null 2>&1
                echo "程序内存占用超过${MAX_MEM_MB}MB限制，已强制终止" >> "$OUTPUT_TMP"
                exit 1
            fi
            sleep 0.1  # 每0.1秒检查一次
        done
    ) &
    MONITOR_PID=$!

    # 执行程序并记录时间和内存信息
    /usr/bin/time -f "%M" "$PROGRAM" "$file" >> "$OUTPUT_TMP" 2>> "$TIME_TMP"
    PROGRAM_EXIT_CODE=$?

    # 停止监控进程
    kill -9 "$MONITOR_PID" >/dev/null 2>&1
    wait "$MONITOR_PID" 2>/dev/null

    # 计算运行时间
    END_TIME=$(date +%s%N)
    DURATION=$(( (END_TIME - START_TIME) / 1000000 ))  # 转换为毫秒

    # 处理内存占用信息
    if [ $PROGRAM_EXIT_CODE -eq 0 ] && [ -s "$TIME_TMP" ]; then
        MEM_KB=$(cat "$TIME_TMP")
        MEM_MB=$(echo "scale=2; $MEM_KB / 1024" | bc)
        MEM_INFO="${MEM_MB}MB"
    else
        if grep -q "已强制终止" "$OUTPUT_TMP"; then
            MEM_INFO="超过${MAX_MEM_MB}MB（被终止）"
        else
            MEM_INFO="获取失败"
        fi
    fi

    # 将临时输出写入日志
    cat "$OUTPUT_TMP" >> "$LOG_FILE"
    echo "运行时间: $DURATION 毫秒" >> "$LOG_FILE"
    echo "最大内存占用: $MEM_INFO" >> "$LOG_FILE"
    echo "----------------------------------------" >> "$LOG_FILE"
    echo "" >> "$LOG_FILE"

    # 清理临时文件
    rm -f "$OUTPUT_TMP" "$TIME_TMP"

    # 控制台输出结果
    if [ $PROGRAM_EXIT_CODE -eq 0 ]; then
        echo "完成处理，耗时: $DURATION 毫秒，内存占用: $MEM_INFO"
    else
        echo "处理中断，耗时: $DURATION 毫秒，内存占用: $MEM_INFO"
    fi

    FILE_NUM=$((FILE_NUM + 1))
done

echo "========================================" >> "$LOG_FILE"
echo "基准测试完成于: $(date)" >> "$LOG_FILE"
echo "所有结果已保存到: $LOG_FILE"

exit 0
