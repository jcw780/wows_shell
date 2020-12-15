rem experimental multithreaded and vectorized version
emcc --bind -o build/shellWasmTV.js shellWasm.cpp --std=c++17 -O3 -s ASSERTIONS=1 -s MODULARIZE=1 -s EXPORT_NAME="ShellWasmTV"^
 -s ENVIRONMENT="web,worker" -s MALLOC=dlmalloc -s ALLOW_MEMORY_GROWTH=1 -s FILESYSTEM=0 -s SINGLE_FILE=0 --closure 1 -flto^
 -msimd128 -v -Wall -Wextra -Wpedantic -pthread