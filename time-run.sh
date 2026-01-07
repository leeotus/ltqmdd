#!/bin/bash

files=("$1"/*)

# 检查传入的目录是否存在
if [ -d "$1" ]; then
    # 遍历目录下的所有文件
    for file in "${files[@]}"; do
        # 仅处理普通文件
        if [ -f "$file" ]; then
            echo "=== current file: $file ==="
            # 使用 /usr/bin/time 统计（注意：不是内置的time）
            # -v 显示详细资源使用（包含内存）
            /usr/bin/time -v ./build/apps/main "$file"
            echo "------------------------"
        fi
    done
else
    echo "错误：传入的参数不是有效的目录"
fi
