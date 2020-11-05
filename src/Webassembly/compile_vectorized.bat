rem vectorized/unthreaded version
emcc --bind -o shellWasmV.wasm.js shellWasm.cpp --std=c++17 -O3 -s ASSERTIONS=1 -s MODULARIZE=1 -s EXPORT_NAME="ShellWasmV"^
 -s ENVIRONMENT="web" -s MALLOC=emmalloc -s ALLOW_MEMORY_GROWTH=1 -s FILESYSTEM=0 -s SINGLE_FILE=1 --closure 1 -flto -msimd128 -v 