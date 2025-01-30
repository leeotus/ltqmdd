# LTQMDDs : Quantum Multiple-Valued Decision Diagrams with Linear Transformations

LTQMDDs (Linearly Transformed Quantum Multiple-Valued Decision Diagrams) are a novel canonical representation in quantum computing, integrating linear transformations into QMDDs to obtain more compact forms of quantum functions, with semantics based on the rearrangement of matrix entries and involving level exchange procedures and linear sifting algorithms for optimization.



## 1.Environment

This repository have been compiled and executed in Ubuntu22.04LTS, 64bit with a machine equipped with Intel Xeon Platinum 8255C 2.50GHz CPU * 8 and 32GB memory.

- gcc/g++ 8+

- cmake > 3.19



## 2.Getting Started

For compiling:

```c++
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

then, you can run the project in the following way:

```shell
sudo chmod +x ./run.sh
./run.sh <path_to_file> 
```

All the files involved in the experiment are saved in this directory: **"circuits/experiments"**

Take "alu4_201.real" file (from [RevLib](http://revlib.org/)) as an example:

```shell
./run.sh ./circuits/experiments/revLib/alu4_201.real
```

Or you can run all files which are saven in the same directory in one command:

```shell
./run.sh ./circuits/experiments/revLib
```

You can run other circuit files as you like. Currently available files are:

- `Real` (e.g. from [RevLib](http://revlib.org))
- `OpenQASM` (e.g. used by [Qiskit](https://github.com/Qiskit/qiskit))
- `TFC` (e.g. from [Reversible Logic Synthesis Benchmarks Page](http://webhome.cs.uvic.ca/~dmaslov/mach-read.html))
- `QC` (e.g. from [Feynman](https://github.com/meamy/feynman))