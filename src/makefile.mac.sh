#!/bin/bash

# 设置源文件和目标文件名
SOURCE="waster_lite_t2.cpp"
TARGET="waster_lite_b2.exe"

# 编译器检测：优先 clang++，其次 g++
if command -v clang++ &>/dev/null; then
    COMPILER="clang++"
elif command -v g++ &>/dev/null; then
    COMPILER="g++"
else
    echo "错误：未找到可用的 C++ 编译器（clang++ 或 g++）" >&2
    exit 1
fi

echo "使用编译器: $COMPILER"

# 执行编译（-std=c++23, -Ofast, -mcmodel=large）
$COMPILER -std=c++23 -Ofast -mcmodel=large "$SOURCE" -o "$TARGET"

# 检查编译是否成功
if [ $? -eq 0 ]; then
    echo "编译成功: $TARGET"
else
    echo "编译失败" >&2
    exit 1
fi