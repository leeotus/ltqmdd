#!/bin/bash

files=("$1"/*)

if [ -d "$1" ]; then
    for file in "${files[@]}"; do
      if [ -f "$file" ]; then
        echo "$file:"
        ./build/apps/main "$file"
      fi
    done
fi
