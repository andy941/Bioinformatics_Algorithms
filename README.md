C++ implementation of the exercises and examples from the book: 

***"Bioinformatics Algorithms - design and implementation in python"*** by Miguel Rocha & Pedro G. Ferreira.

I wanted to know more about algorithms to sharpen up my bioinformatics skills. It seemed fitting to follow the example and exercises in C++ having just finished a [comprehensive introduction in the language](https://github.com/andy941/Cpp_principles_exercises#cpp_principles_exercises).

---

#### Build project binaries:
``` bash
mkdir build/
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cd ..
cmake --build build
cmake --install build
```

#### Additional setup:
Eigen matrix library is all headers, I added it as a git submodule. Just make sure is downloaded properly.
``` bash
git submodule update --init --recursive
```
