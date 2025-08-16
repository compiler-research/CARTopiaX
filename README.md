# CARTopiaX

This repository provides an **agent-based simulation** of tumour-derived organoids and their interaction with **CAR T-cell therapy**.  
Developed as part of Google Summer of Code 2025, the project is released under the **Apache License 2.0**.

The simulation integrates computational modeling and biological insights to explore tumour–immune dynamics and assess treatment outcomes under various scenarios.

---

## Table of Contents

1. [Project Overview](#project-overview)
2. [Dependencies](#dependencies)
3. [Installation](#installation)
4. [Building the Simulation](#building-the-simulation)
5. [Running the Simulation](#running-the-simulation)
6. [Acknowledgments](#acknowledgments)
7. [License](#license)

---

## Project Overview

This project implements an **agent-based model** to simulate the behavior of *tumour-derived organoids* (lab-grown tumour models) and their response to **CAR T-cell therapy**.

With this simulation, researchers can:
- Recreate in vitro conditions for tumour growth.
- Introduce CAR T-cells and analyze their effectiveness in solid tumor microenvironments.
- Explore different treatment strategies and parameter variations.
- Evaluate outcomes such as tumour reduction, elimination, or relapse risk.

By adjusting biological and therapeutic parameters, the model enables **in silico experimentation** to complement laboratory research.

---

## Dependencies

- [BioDynaMo](https://biodynamo.org/) (tested with version 1.05.132)
- CMake ≥ 3.13
- GCC or Clang with C++17 support
- GoogleTest (for unit testing)

**Note:** Ensure BioDynaMo is installed and sourced before running the simulation.

---

## Installation

Clone the repository:
```bash
git clone https://github.com/compiler-research/CARTopiaX.git
cd cart_tumor
```

---

## Building the Simulation

**Option 1:**  
Use BioDynaMo’s build system:
```bash
biodynamo build
```

**Option 2:**  
Manual build:
```bash
mkdir build && cd build
cmake ..
make -j <number_of_processes>
```

---

## Running the Simulation

After building, run the simulation using one of the following methods:

**Option 1:**  
With BioDynaMo:
```bash
biodynamo run
```

**Option 2:**  
Directly from the build directory:
```bash
./build/cart_tumor
```

---

## Acknowledgments

This project builds upon the BioDynaMo simulation framework.

> Lukas Breitwieser, Ahmad Hesam, Jean de Montigny, Vasileios Vavourakis, Alexandros Iosif, Jack Jennings, Marcus Kaiser, Marco Manca, Alberto Di Meglio, Zaid Al-Ars, Fons Rademakers, Onur Mutlu, Roman Bauer.  
> *BioDynaMo: a modular platform for high-performance agent-based simulation*.  
> Bioinformatics, Volume 38, Issue 2, January 2022, Pages 453–460.  
> [https://doi.org/10.1093/bioinformatics/btab649](https://doi.org/10.1093/bioinformatics/btab649)

Some of the mathematical models and solver implementations are based on the research of  
Luciana Melina Luque and collaborators, as described in:

> Luque, L.M., Carlevaro, C.M., Rodriguez-Lomba, E. et al.  
> *In silico study of heterogeneous tumour-derived organoid response to CAR T-cell therapy*.  
> Scientific Reports 14, 12307 (2024).  
> [https://doi.org/10.1038/s41598-024-63125-5](https://doi.org/10.1038/s41598-024-63125-5)

---

## License

This project is licensed under the Apache License 2.0. See the [LICENSE](LICENSE) file for details.