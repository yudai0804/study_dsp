# デジタル信号処理(DSP)勉強用リポジトリ

# cmakeメモ
cmakeを使ってMakefileを生成
```
cd path/to/target/directory
```
cmakeコマンドを使って、makefileを生成する  
そのとき、-SオプションでSourceディレクトリを、-Bオプションでbuildディレクトリを生成する  
```
cmake -S . -B build
```
以下のコマンドでビルドする。  
使うコマンドは好きなほうで。
```
cd build
make
```
or
```
cmake --build build
```
これで実行ファイルが生成される。