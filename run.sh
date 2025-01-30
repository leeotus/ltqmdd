#!/bin/bash

if [ $# -ne 1 ]; then
    echo "Usage: $0 <path_to_file> or <path_to_folder>"
    exit 1
fi

files=("$1"/*)

if [ -d "$1" ]; then
    for file in "${files[@]}"; do
        if [ -f "$file" ]; then
            echo "$file:"
            ./build/apps/orgnl "$file" &
            wait
            ./build/apps/mixed "$file" &
            wait
            ./build/apps/upper "$file" &
            wait
            ./build/apps/lower "$file" &
            wait
        fi
    done
    exit 1
fi

./build/apps/orgnl "$1" &
./build/apps/mixed "$1" &
./build/apps/upper "$1" &
./build/apps/lower "$1" &

wait
