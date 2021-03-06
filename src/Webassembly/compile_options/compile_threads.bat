rem experimental multithreaded version
emcc --bind -o ../build/shellWasmT.js ../shellWasm.cpp --std=c++17 -O3 -s ASSERTIONS=1 -s MODULARIZE=1 -s EXPORT_NAME="ShellWasmT"^
 -s ENVIRONMENT="web,worker" -s MALLOC=dlmalloc -s ALLOW_MEMORY_GROWTH=1 -s FILESYSTEM=0 -s SINGLE_FILE=0 --closure 1 -flto^
 -v -Wall -Wextra -Wpedantic -pthread -s EXPORT_ES6=1 -s USE_ES6_IMPORT_META=0