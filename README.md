# CARTopiaX

<table>
  <tr>
    <td>
      <div style="background-color:white; padding:10px; display:inline-block;">
        <img src="./assets/images/CARTopiaXLogo_4x.png" width="200" alt="CARTopiaX logo">
      </div>
    </td>
    <td>
      <p>
        This repository provides an <strong>agent-based simulation</strong> of tumor-derived organoids and their interaction with <strong>CAR T-cell therapy</strong>.<br>
        Developed as part of Google Summer of Code 2025 at the CERN HEP Software Foundation (HSF), the project is released under the <strong>Apache License 2.0</strong>.<br>
        <br>
        The simulation integrates computational modeling and biological insights to explore tumor–immune dynamics and assess treatment outcomes under various scenarios.
      </p>
    </td>
  </tr>
</table>


---

## Table of Contents

1. [Project Overview](#project-overview)
2. [Project Structure](#project-structure)
3. [Model Replication](#model-replication)
4. [Dependencies](#dependencies)
5. [Installation](#installation)
6. [Development Environment (Container)](#development-environment-container)
7. [Building the Simulation](#building-the-simulation)
8. [Input Parameters](#input-parameters)
9. [Running the Simulation](#running-the-simulation)
10. [Visualizing Results](#visualizing-results)
11. [Model Calibration (Bayesian Optimization)](#model-calibration-bayesian-optimization)
12. [Acknowledgments](#acknowledgments)
13. [License](#license)
14. [Author Contact Information](#author-contact-information)


---

## Project Overview

**CAR T-cell therapy** is a form of cancer immunotherapy that engineers a patient’s T cells to recognize and destroy malignant cells. While this approach has achieved remarkable success in treating blood cancers, it faces significant challenges in solid tumors due to the complexity and heterogeneity of their microenvironments.

**CARTopiaX** is an **agent-based model** designed to simulate the behavior of *tumour-derived organoids* which are lab-grown models that mimic real solid tumor environments and their response to CAR T-cell therapy.  
The project aims to bridge the gap between laboratory experiments and computational biology by developing a high-fidelity **in silico digital twin** of tumor-derived organoids for studying solid tumor dynamics and immunotherapy outcomes.

Built on **BioDynaMo**, a high-performance, open-source platform for large-scale biological modeling, CARTopiaX enables researchers to:

- Recreate realistic *in vitro* conditions for tumor growth.  
- Introduce CAR T-cells and evaluate their efficacy in heterogeneous, solid tumor microenvironments.  
- Explore different therapeutic strategies and parameter variations.  
- Assess treatment outcomes such as tumor reduction, elimination, or relapse risk.  

CARTopiaX implements the mathematical framework described in the *Nature* publication *“In silico study of heterogeneous tumour-derived organoid response to CAR T-cell therapy,”* successfully replicating its key results and extending them through improved scalability and performance.

### Project Highlights

- **Performance:** Simulations run more than twice as fast as previous existing model, enabling rapid scenario exploration and hypothesis validation.  
- **Software Quality:** Developed in **C++** following robust software engineering practices, ensuring high-quality, maintainable, and efficient code.  
- **Architecture:** Designed to be scalable, modular, and extensible, fostering collaboration, customization, and continuous evolution within an open-source ecosystem.  

Together, these features make **CARTopiaX** a powerful computational tool for investigating CAR T-cell dynamics in solid tumors by accelerating scientific discovery, guiding experimental design, and reducing the cost and time associated with wet-lab research.

---

## Project Structure

The diagram below illustrates the **agents and main mechanics** interacting within CARTopiaX:


<!-- ![Model Diagram](./assets/images/model_outline.png) -->
```mermaid
---
config:
  flowchart:
    curve: basis
---
flowchart LR
    CART("CAR-T Cell")
    IMMUNO(("Immunostimulatory<br>factor"))
    OXY(("Oxygen"))
    GLU(("Glucose"))
    TUMOR("Tumor Cell")
    ONCO(("Oncoprotein<br>Level"))
    NOTE["Represents<br>the antigen<br>expression of the<br>tumor cell.<br>Cells with higher<br>oncoprotein levels are<br>more proliferative but<br>also more likely to be<br>killed by CAR-T cells."]

    CART -->|"moved by<br>chemotaxis<br>towards"| IMMUNO
    IMMUNO <-->|"secretes"| TUMOR
    CART -->|"Induces apoptosis to"| TUMOR
    OXY -->|"lack of it<br>induces<br>necrosis"| TUMOR
    OXY <-->|"consumes"| TUMOR
    CART -->|"consumes"| OXY
    GLU <-->|"consumes"| TUMOR
    GLU -->|"lack of it<br>induces<br>apoptosis"| TUMOR
    CART -->|"consumes"| GLU
    TUMOR -->|"expresses"| ONCO
    ONCO -.- NOTE

    linkStyle 1,4,6 marker-end:none
    classDef cell font-size:22px
    class CART,TUMOR cell
```



The project is organized into the following components:

### Core Simulation Files (`src/`)

- **[`tumor_cell.h`](src/agents/tumor_cell.h) / [`tumor_cell.cc`](src/agents/tumor_cell.cc)**: Defines tumor cells with four states (alive, necrotic swelling, necrotic lysed, apoptotic), four types based on oncoprotein levels, volume dynamics using exponential relaxation, oxygen/immunostimulatory factor exchange, and the [`StateControlGrowProliferate`](src/agents/tumor_cell.h) behavior for cancer growth and state transitions.

- **[`cart_cell.h`](src/agents/cart_cell.h) / [`cart_cell.cc`](src/agents/cart_cell.cc)**: Implements CAR-T cells with alive/apoptotic states, tumor cell attachment/detachment mechanisms, stochastic killing attempts, chemotaxis toward immunostimulatory factors, finite lifespan, and the [`StateControlCart`](src/agents/cart_cell.h) behavior.

- **[`cart_tumor.h`](src/cart_tumor.h) / [`cart_tumor.cc`](src/cart_tumor.cc)**: Contains the main [`Simulate`](src/cart_tumor.cc) function that configures BioDynaMo, creates diffusion grids for oxygen and immunostimulatory factors, initializes the tumor sphere, sets up custom mechanical forces, and schedules treatment/output operations.

- **[`forces_tumor_cart.h`](src/forces/forces_tumor_cart.h) / [`forces_tumor_cart.cc`](src/forces/forces_tumor_cart.cc)**: Implements the [`InteractionVelocity`](src/forces/forces_tumor_cart.h) class with velocity-dependent repulsion and adhesion forces between cells, with differentiated force coefficients for tumor-tumor, CAR-CAR, and tumor-CAR interactions.

- **[`diffusion_thomas_algorithm.h`](src/diffusion/diffusion_thomas_algorithm.h) / [`diffusion_thomas_algorithm.cc`](src/diffusion/diffusion_thomas_algorithm.cc)**: Solves the 3D diffusion equation (∂t u = ∇D∇u - μu) for chemical substances (oxygen and immunostimulatory factors) using the Thomas algorithm and Alternating Direction Implicit (ADI) method for tridiagonal systems in each spatial direction. Supports Dirichlet boundary conditions and integrates cellular consumption/secretion (∂ρ/∂t = ∇·(D∇ρ) − λ·ρ + Σ[(V_k/V_voxel)·(S_k·(ρ*_k − ρ) − (S_k + U_k)·ρ)]) through [`ComputeConsumptionsSecretions`](src/diffusion/diffusion_thomas_algorithm.cc), which updates concentrations based on cell-specific uptake and secretion rates.

- **[`hyperparams.h`](src/params/hyperparams.h) / [`hyperparams.cc`](src/params/hyperparams.cc)**: Contains the [`SimParam`](src/params/hyperparams.h) class with default biological/simulation parameters and treatment schedule definition. Provides [`LoadParams`](src/params/hyperparams.cc) to read configuration from [`params.json`](params.json) and [`PrintParams`](src/params/hyperparams.cc) to display current parameter values.

- **[`utils_aux.h`](src/utils/utils_aux.h) / [`utils_aux.cc`](src/utils/utils_aux.cc)**: Provides utility functions ([`SamplePositiveGaussian`](src/utils/utils_aux.cc), [`CreateSphereOfTumorCells`](src/utils/utils_aux.cc), [`AnalyzeTumor`](src/utils/utils_aux.cc), [`GenerateRandomDirection`](src/utils/utils_aux.cc)) and operations ([`SpawnCart`](src/utils/utils_aux.cc) for dosage administration, [`OutputSummary`](src/utils/utils_aux.cc) for CSV data export to `/output/final_data.csv`).

- **[`substance_interactor.h`](src/interfaces/substance_interactor.h)**: Defines the [`ISubstanceInteractor`](src/interfaces/substance_interactor.h) virtual interface, which must be implemented by any agent that consumes or secretes diffusing substances (oxygen and immunostimulatory factor). Both [`TumorCell`](src/agents/tumor_cell.h) and [`CarTCell`](src/agents/cart_cell.h) implement this interface via the [`ConsumeSecreteSubstance`](src/interfaces/substance_interactor.h) method, allowing the diffusion grid's [`ComputeConsumptionsSecretions`](src/diffusion/diffusion_thomas_algorithm.cc) to dispatch to the correct cell-type logic.

### Configuration Files

- **[`params.json`](params.json)**: JSON-formatted parameter file for configuring simulation runs without recompilation. In [`params4SphericalModel.json`](params4SphericalModel.json) there is a configuration example for an spherical tumor and in [`params4CylindricalModel.json`](params4CylindricalModel.json) the configuration to simulate a cylindrical one.

- **[`bdm.toml`](bdm.toml)**: BioDynaMo-specific configuration for visualization settings.

- **[`CMakeLists.txt`](CMakeLists.txt)**: CMake build configuration that integrates with BioDynaMo and sets up the C++17 project structure. Modification is typically not required.

### Analysis

- **[`CARTopiaX_Simulation_Analysis.ipynb`](CARTopiaX_Simulation_Analysis.ipynb)**: Jupyter notebook for post-processing simulation results, generating plots, and statistical analysis.

### Model Calibration (`abm_calibration/`)

- **[`optimize.py`](abm_calibration/optimize.py)**: Bayesian optimization workflow to calibrate model parameters against target data. See [Model Calibration](#model-calibration-bayesian-optimization).

- **[`requirements.txt`](abm_calibration/requirements.txt)**: Python dependencies for the calibration workflow.


---

## Model Replication

All plots in this section compare CARTopiaX with the model results from the *Nature* paper over five runs, demonstrating **successful replication** and reproducing the same **biological findings**. Even though the graphs do not always overlap, this is due to substantial known differences in their modeling approaches and stochastic nature. What is important is that the overall behaviors and key biological dynamics are accurately reproduced, as researchers primarily focus on these trends and peak responses when designing treatments.

### Replication Example: Tumor with no CAR-T treatment
---
The following example reproduces the 30-day evolution of a 150 µm radius tumor simulation with no CAR T-cell treatment applied. The plot below shows the total number of cancer cells and the tumor radius over time:

![No treatment cells and radius](./assets/images/no_doses_num_cells_and_tumor_radius.png)

#### Tumor Heterogeneity:

CARTopiaX models a **heterogeneous tumor**, where cancer cells differ in their *oncoprotein expression levels*, representing varying proliferative capacities.  
Although oncoprotein levels are continuous, cells are grouped into four discrete categories:

| Type | Oncoprotein Level     | Aggressiveness       |
|:----:|:---------------------:|:--------------------:|
| 1    | 1.5 < Oncoprotein ≤ 2 | Most proliferative   |
| 2    | 1 < Oncoprotein ≤ 1.5 | Very proliferative   |
| 3    | 0.5 < Oncoprotein ≤ 1 | Less proliferative   |
| 4    | 0 ≤ Oncoprotein ≤ 0.5 | Least proliferative  |

Tumor evolution visualized in ParaView:

![No treatment ParaView](./assets/images/no_dose_tumor_evolution_paraview.png)

Type 1 cells are highly proliferative, causing their proportion in the tumor to increase, while Type 3–4 cells divide less frequently and, struggling in a oxygen-limited, resource-competitive environment, gradually become less common:

![No treatment Type of Cells](./assets/images/no_doses_type_of_cells_absolute_numbers.png)

As cell quantity rises, **oxygen levels drop** due to increased consumption, while the average **oncoprotein level gradually increases** because more proliferative cells pass their higher oncoprotein levels to their progeny over time.

![No treatment Oxygen and Oncoprotein](./assets/images/no_doses_oncoprotein_and_oxygen.png)


### Replication Example: One CAR T-cell dose of Scale 1:1
---
The following example models the same 150 µm-radius tumor, but with a single CAR T-cell dose, equal in number to the tumor cells administered, on day 0. The CAR T-cell population then declines stochastically due to apoptosis, reaching minimal levels around day 10. The plot below shows the total tumor and CAR T-cell populations over time:

![1Dose Tumor and CAR T-cell amount](./assets/images/dose_scale1_day0_num_cells.png)

All cells consume oxygen, and cancer cells can die from **necrosis** in its absence. In addition, CAR T-cells tend to **eliminate higher oncoprotein-expressing cells more effectively**, leading to the following dynamics:

<u>Before day 10:</u> CAR T-cells are actively killing tumor cells and **Type 1 and 2 cell populations decrease** as more proliferative cancer cells are targeted.

<u>After day 10:</u> When CAR T-cells die from apoptosis, **Type 1 and 2 cells increase** their proportion in the tumor at the expense of Type 3 and 4 cells, as high-oncoprotein expressers divide faster.

![1dose Percentage of Each Cell Type](./assets/images/dose_scale1_day0_type_of_cells_percentage_populations.png)

<u>Before day 10:</u> In the beginning, **oxygen** levels decrease with CAR T-cell arrival. Then they **gradually rise** as both CAR T-cells and tumor cells die, reducing overall consumption. On the other hand, the **average oncoprotein level drops** rapidly due to CAR T-cell preferential elimination of the most aggressive cancer cells.

<u>After day 10:</u> CAR T-cells are no longer present and therefore tumor resumes growth, making **oxygen levels decline**. In addition, **oncoprotein levels rise** as highly proliferative cells dominate once again.

![1dose Oncoprotein and Oxygen Levels](./assets/images/dose_scale1_day0_oncoprotein_and_oxygen.png)

During the simulation, **dead and resistant cells** accumulate around the tumor core, forming a ***shield-like barrier*** that impedes CAR T-cell infiltration and reduces treatment effectiveness.

![Cell Shield Formation](./assets/images/shield.png)

A 3D visualization of the sliced tumor and CAR T-cell treatment used in this example, rendered in *ParaView*, is available here:  
[Watch the simulation video](https://youtu.be/7V8n627Nmzc)

### Replication Example: Less is better, increasing cellular dosage does not always increase efficacy

---
We consider again the 150 µm-radius tumor and now **two CAR T-cell doses of scale 1:1** (equal in number to the initial tumor cells) are delivered on days 0 and 8. This way much less tumor cells are present at day 30 in comparison to applying just one dose:

![2Doses of Scale 1:1](./assets/images/dose_scale1_day0_dose_scale1_day8_num_cells.png)

Delivering more CAR T-cells has helped in controlling the tumor, however this is not always the case. The next plot shows the same initial tumor but with two CAR T-cell doses containing a quantity twice the initial tumor cell count, on days 0 and 5. By day 30, the number of **tumor cells is roughly the same** as before, despite using **twice the amount of CAR T-cells**.

![2Doses of Scale 2:1](./assets/images/dose_scale2_day0_dose_scale2_day5_num_cells.png)


Increasing CAR T-cell dosage does not necessarily improve tumor killing and can increase *toxicity*. The model suggests two doses at a 1:1 CAR T-to-cancer cell ratio, balancing effectiveness and safety, minimizing inactive *free* CAR T-cells.

### Performance
---

CARTopiaX successfully replicates the findings described in the *Nature* publication, achieving a **2× speed improvement** compared to the previous implementation.  
This performance gain enables faster scenario exploration and larger-scale simulations.

![Execution Time Comparison](./assets/images/execution_times.png)

### Scripts Proving Model Replication
---

This [repository](https://github.com/salva24/CARTopiaX_scripts_proving_model_replication) contains the scripts used to validate and compare simulation results between CARTopiaX and its foundational model, CART-ABM.
It can be employed to verify that the reported improvements are objectively supported by the data.

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
cd CARTopiaX

```

---

## Development Environment (Container)

A prebuilt BioDynaMo — with its bundled ROOT — is published as a
[ci-workflows](https://github.com/compiler-research/ci-workflows) recipe cell,
so you do not have to build it yourself. `bin/start` drops you into a persistent
container with that BioDynaMo in place, your checkout mounted read-write, and
[Claude Code](https://claude.com/claude-code) installed.

This is the fastest way for a new contributor to get a working environment.

### Prerequisites

- **Docker**, running.
- **git** and **Python 3** on the host.
- **~3 GB free** in the Docker VM (~250 MB download, ~800 MB unpacked, plus build space).

`nektos/act` is *not* needed. It is only required to run a whole CI job locally
with `bin/repro <row-name>`, which asks `act` to expand the matrix.

> On Apple Silicon the container is x86_64 and runs translated. That is fine for
> developing, but do not trust it for timing measurements or for diagnosing
> toolchain-level failures.

### Getting in

```bash
git clone https://github.com/compiler-research/ci-workflows ~/sources/ci-workflows
cd ~/sources/ci-workflows && ./bin/start
```

Pick CARTopiaX from the list. It clones the project, downloads the BioDynaMo
build CI uses, installs its host dependencies and opens a shell. Already have a
checkout? Run `~/sources/ci-workflows/bin/start` from inside it and the menu is
skipped.

The first entry downloads and unpacks BioDynaMo; later entries take seconds.
`claude` is installed for you — log in once per container. A host that has never
run Claude Code needs nothing extra; if you *do* already run it, anything you
put in `~/.cache/ci-workflows/devshell-cache/ai/skills` is what it sees inside.

### Building, inside the container

```bash
source "$DEVSHELL_INSTALL/bin/thisbdm.sh"

cd /patches
cmake . -B build -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++
cmake --build build -j"$(nproc)"
ctest --test-dir build --output-on-failure
```

The compiler pin is a workaround rather than a preference, so do not drop it:
the devshell exports `CC=clang`, while BioDynaMo replaces the compiler with
MPI's wrapper (which wraps gcc) *after* CMake has detected clang. The OpenMP
flags then disagree and g++ rejects `-fopenmp=libomp`.

### Modifying BioDynaMo itself

The recipe ships BioDynaMo's source next to the install, at `$DEVSHELL_SRC` —
a shallow checkout at the tag CI pins, already configured at `$DEVSHELL_BUILD`,
so `make -C $DEVSHELL_BUILD` iterates directly.

To build CARTopiaX against your modified BioDynaMo, install it to a prefix of
your own and source *that* instead of `$DEVSHELL_INSTALL`:

```bash
cmake --install "$DEVSHELL_BUILD" --prefix ~/bdm-dev
source ~/bdm-dev/biodynamo-v1.05/bin/thisbdm.sh
```

The install is self-contained — it carries its own ROOT — and `thisbdm.sh`
derives `BDMSYS` from its own location, so it works from anywhere. Note the
extra `biodynamo-v<version>` level: BioDynaMo installs one directory deeper than
the prefix you give it, whereas `$DEVSHELL_INSTALL` is that inner directory
already, flattened by the recipe. Leaving `$DEVSHELL_INSTALL` untouched keeps
the pristine CI environment one `source` away.

### What persists

| layer | holds | survives |
|---|---|---|
| container `devshell-<cell>` | apt packages, Claude login, your `$HOME` | exiting the shell |
| host cache `~/.cache/ci-workflows/devshell-cache` | BioDynaMo install and source, ccache, Claude skills/settings/memory | removing the container |
| your checkout | **your work** | everything — it is your working copy |

Your work is never inside the container: `/patches` *is* your checkout, so
commit, push and open pull requests from the host as usual. No GitHub
credentials are copied into the container.

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

## Input Parameters

All hyperparameters are modifiable by giving them a value in the `params.json` without the need to recompile the model after changing them. When a value is not specified the default values are taken. Some hyperparameters are derived from others by default but they can still be set to arbitrary values if wished. Don't worry if you modify a hyperparameter in `params.json` but not its derived values; by default, derived values are computed at the start of the simulation based on the specified ones.

### General Simulation Parameters

These are basic parameters that are commonly changed when designing a treatment with CARTopiaX.


| Parameter | Default Value | Units | Description |
|-----------|---------------|-------|-------------|
| `seed` | 42 | - | Seed for random number generation to ensure reproducibility |
| `output_performance_statistics` | false | - | Enable/disable performance statistics output |
| `total_minutes_to_simulate` | 43200 | minutes | Total simulation time (default: 30 days) |
| `tumor_shape` | `"sphere"` | - | Shape of the initial tumor: `"sphere"` or `"cylinder"` |
| `initial_spherical_tumor_radius` | 150 | μm | Initial radius of the spherical tumor (only used if `tumor_shape` is `"sphere"`) |
| `cylindrical_tumor_radius` | 3000 | μm | Initial radius of the cylindrical tumor (only used if `tumor_shape` is `"cylinder"`) |
| `cylindrical_tumor_height` | 100 | μm | Initial height of the cylindrical tumor (only used if `tumor_shape` is `"cylinder"`) |
| `initial_number_of_cylindrical_tumor_cells` | 2800 | cells | Initial number of tumor cells in the cylindrical tumor (only used if `tumor_shape` is `"cylinder"`) |
| `bounded_space_length` | 1000 | μm | Length of the cubic simulation domain |
| `bound_space_toplogy` | `"torus"` | - | Topology of the domain boundaries for the agents: `"torus"`, `"open"` or `"closed"` (see BioDynaMo documentation) |
| `treatment` | `{0: 3957, 8: 3957}` | day:cells | Map of treatment days to CAR-T cell counts. Key = day, Value = number of cells |

### Specific Parameters: 
These are all the other hyperparameters that are not adviced to be modified unless having a deeper understanding of the model and code.


#### Domain Boundary Parameters

These parameters restrict the region where cells can be located, e.g. to model physical barriers in cylindrical tumors. Coordinates are expressed in the simulation coordinate system, where the domain extends from `-bounded_space_length/2` to `+bounded_space_length/2`. By default they do not impose any restriction.

| Parameter | Default Value | Units | Description |
|-----------|---------------|-------|-------------|
| `bounded_space_min_allowed_z` | `-bounded_space_length/2` (-500) | μm | Minimum allowed z coordinate for cells |
| `bounded_space_max_allowed_z` | `bounded_space_length/2` (500) | μm | Maximum allowed z coordinate for cells |
| `bounded_space_max_allowed_radius` | `bounded_space_length` (1000) | μm | Maximum allowed distance of cells from the center of the tumor (from the center point for spherical tumors and from the axis for cylindrical ones) |


#### Time Step Parameters

| Parameter | Default Value | Units | Description |
|-----------|---------------|-------|-------------|
| `dt_substances` | 0.01 | minutes | Time step for diffusion and substance exchange |
| `dt_mechanics` | 0.1 | minutes | Time step for mechanical forces between cells |
| `dt_cycle` | 6 | minutes | Time step for cell cycle progression |
| `dt_step` | 0.1 | minutes | General simulation time step (same as `dt_mechanics`) |

#### Output Parameters

| Parameter | Default Value | Units | Description |
|-----------|---------------|-------|-------------|
| `output_csv_interval` | 7200 | steps | Interval for writing simulation data to CSV (default: 12 hours) |
| `output_information_dependent_on_radius` | false | - | Outputs an additional CSV (`output/data_dependent_on_radius_tumor.csv`) with information aggregated by distance from the center of the tumor (from the center point for spherical tumors and from the axis for cylindrical ones) |
| `max_radius_analysis_csv_dependent_on_radius` | 3000 | μm | Maximum radius considered in the radius-dependent CSV |
| `num_radius_intervals` | 10 | - | Number of equal radius intervals, from 0 to `max_radius_analysis_csv_dependent_on_radius`, in which the radius-dependent information is aggregated |

#### Apoptosis Parameters

| Parameter | Default Value | Units | Description |
|-----------|---------------|-------|-------------|
| `volume_relaxation_rate_cytoplasm_apoptotic_cells` | 0.0166667 | min⁻¹ | Cytoplasm volume relaxation rate for apoptotic cells |
| `volume_relaxation_rate_nucleus_apoptotic_cells` | 0.00583333 | min⁻¹ | Nucleus volume relaxation rate for apoptotic cells |
| `volume_relaxation_rate_fluid_apoptotic_cells` | 0.0 | min⁻¹ | Fluid volume relaxation rate for apoptotic cells |
| `time_apoptosis` | 516 | minutes | Time until an apoptotic cell is removed from simulation |
| `reduction_consumption_dead_cells` | 0.1 | - | Reduction factor for oxygen consumption in dead cells |

#### Chemical Diffusion Parameters

| Parameter | Default Value | Units | Description |
|-----------|---------------|-------|-------------|
| `resolution_grid_substances` | 50 | voxels/axis | Number of voxels per axis for diffusion grids |
| `min_initial_z_substances` | `-bounded_space_length/2` (-500) | μm | Minimum height at which the substance grids are initialized to a value different from 0 |
| `max_initial_z_substances` | `bounded_space_length/2` (500) | μm | Maximum height at which the substance grids are initialized to a value different from 0 |

<u>Oxygen Parameters</u>

| Parameter | Default Value | Units | Description |
|-----------|---------------|-------|-------------|
| `diffusion_coefficient_oxygen` | 100000 | μm²/min | Diffusion coefficient for oxygen |
| `decay_constant_oxygen` | 0.1 | min⁻¹ | Decay constant (λ) for oxygen |
| `oxygen_reference_level` | 38 | mmHg | Boundary condition value for oxygen concentration |
| `initial_oxygen_level` | 38 | mmHg | Initial oxygen concentration in all voxels |
| `lateral_oxygen_production_min_z` | `-bounded_space_length/2` (-500) | μm | Minimum height at which the Dirichlet oxygen boundary condition is applied. If it equals the domain's minimum z, the floor also produces oxygen |
| `lateral_oxygen_production_max_z` | `bounded_space_length/2` (500) | μm | Maximum height at which the Dirichlet oxygen boundary condition is applied. If it equals the domain's maximum z, the roof also produces oxygen |
| `diffuse_oxygen_on_z_axis` | true | - | Whether oxygen also diffuses along the z axis. If false, it only diffuses in the x and y axes (ideal for a 2D model) |

<u>Immunostimulatory Factor Parameters</u>

| Parameter | Default Value | Units | Description |
|-----------|---------------|-------|-------------|
| `add_immunostimulatory_factor` | true | - | Whether to add the immunostimulatory factor diffusion grid |
| `diffusion_coefficient_immunostimulatory_factor` | 1000 | μm²/min | Diffusion coefficient for immunostimulatory factor |
| `decay_constant_immunostimulatory_factor` | 0.016 | min⁻¹ | Decay constant (λ) for immunostimulatory factor |
| `diffuse_immunostimulatory_factor_on_z_axis` | true | - | Whether the immunostimulatory factor also diffuses along the z axis. If false, it only diffuses in the x and y axes (ideal for a 2D model) |

<u>Glucose Parameters</u>

Glucose is disabled by default. When `add_glucose` is set to `true`, a glucose diffusion grid is added and glucose levels affect tumor cell growth and death (see the glucose parameters in [Tumor Cell Parameters](#tumor-cell-parameters)).

| Parameter | Default Value | Units | Description |
|-----------|---------------|-------|-------------|
| `add_glucose` | false | - | Whether to add the glucose diffusion grid |
| `diffusion_coefficient_glucose` | 7800 | μm²/min | Diffusion coefficient for glucose |
| `decay_constant_glucose` | 0.01 | min⁻¹ | Decay constant (λ) for glucose |
| `initial_glucose_level` | 24.98 | mmol/L | Initial glucose concentration in each voxel |
| `max_radius_glucose_initialization` | `bounded_space_length` (1000) | μm | Radius of the sphere (spherical tumor) or cylinder (cylindrical tumor) inside which glucose is initialized. Outside of it glucose is initialized to 0 |
| `diffuse_glucose_on_z_axis` | true for `"sphere"`, false for `"cylinder"` | - | Whether glucose also diffuses along the z axis. If false, it only diffuses in the x and y axes (ideal for a 2D model) |

#### Mechanical Forces Parameters

| Parameter | Default Value | Units | Description |
|-----------|---------------|-------|-------------|
| `cell_repulsion_between_tumor_tumor` | 10 | - | Repulsion coefficient between tumor cells |
| `cell_repulsion_between_cart_cart` | 50 | - | Repulsion coefficient between CAR-T cells |
| `cell_repulsion_between_cart_tumor` | 50 | - | Repulsion coefficient from CAR-T to tumor cells |
| `cell_repulsion_between_tumor_cart` | 10 | - | Repulsion coefficient from tumor to CAR-T cells |
| `max_relative_adhesion_distance` | 1.25 | - | Maximum relative distance for adhesion (multiplier of cell radius) |
| `cell_adhesion_between_tumor_tumor` | 0.4 | - | Adhesion coefficient between tumor cells |
| `cell_adhesion_between_cart_cart` | 0 | - | Adhesion coefficient between CAR-T cells |
| `cell_adhesion_between_cart_tumor` | 0 | - | Adhesion coefficient from CAR-T to tumor cells |
| `cell_adhesion_between_tumor_cart` | 0 | - | Adhesion coefficient from tumor to CAR-T cells |
| `length_box_mechanics` | 22 | μm | Box length for spatial partitioning in force calculations |
| `dnew` | 0.15 | - | Adams-Bashforth coefficient for current velocity (dt × 1.5) |
| `dold` | -0.05 | - | Adams-Bashforth coefficient for previous velocity (dt × -0.5) |

#### Tumor Cell Parameters

| Parameter | Default Value | Units | Description |
|-----------|---------------|-------|-------------|
| `oncoprotein_mean` | 1 | - | Mean oncoprotein expression level in tumor cells |
| `oncoprotein_standard_deviation` | 0.25 | - | Standard deviation of oncoprotein expression |
| `oncoprotein_limit` | 0.5 | - | Minimum oncoprotein level for CAR-T recognition |
| `oncoprotein_saturation` | 2.0 | - | Maximum oncoprotein level |
| `time_lysis` | 86400 | minutes | Time until lysed necrotic cell removal |
| `basal_death_probability_cancer_cells` | 0 | min⁻¹ | Basal death probability of tumor cells due to natural random causes |
| `default_volume_new_tumor_cell` | 2494 | μm³ | Mean total volume of newly created tumor cell |
| `std_volume_new_tumor_cell` | 0 | μm³ | Standard deviation of the total volume of newly created tumor cells |
| `min_volume_new_tumor_cell` | 2494 | μm³ | Minimum total volume of newly created tumor cells (sampled volumes are clipped) |
| `max_volume_new_tumor_cell` | 2494 | μm³ | Maximum total volume of newly created tumor cells (sampled volumes are clipped) |
| `default_fraction_of_volume_for_nucleus_tumor_cell` |  0.21652 | - | Nuclear volume of newly created tumor cell |
| `default_fraction_fluid_tumor_cell` | 0.75 | - | Fraction of cytoplasmic volume that is fluid |
| `average_time_transformation_random_rate` | 38.6 | hours | Mean cell cycle duration |
| `standard_deviation_transformation_random_rate` | 3.7 | hours | Standard deviation of cell cycle duration |
| `minimum_tumor_cell_target_volume_fraction_for_division` | 0 | - | Minimum fraction of its target volume a tumor cell must reach to be able to divide |
| `adhesion_time` | 60 | minutes | Average time tumor cell remains attached to CAR-T before escaping |

<u>Tumor Cell Oxygen Parameters</u>

| Parameter | Default Value | Units | Description |
|-----------|---------------|-------|-------------|
| `default_oxygen_consumption_tumor_cell` | 10 | 1/min | Baseline oxygen consumption rate of tumor cells |
| `oxygen_saturation_for_proliferation` | 38 | mmHg | Oxygen level for maximum proliferation rate |
| `oxygen_limit_for_proliferation` | 10 | mmHg | Minimum oxygen level for cell proliferation |
| `oxygen_limit_for_necrosis` | 5 | mmHg | Oxygen level below which necrosis begins |
| `oxygen_limit_for_necrosis_maximum` | 2.5 | mmHg | Oxygen level for maximum necrosis probability |
| `maximum_necrosis_lack_of_oxygen_rate` | 0.00277778 | min⁻¹ | Maximum necrosis rate at 0 oxygen (1/360 min⁻¹) |

<u>Tumor Cell Immunostimulatory Factor Parameters</u>

| Parameter | Default Value | Units | Description |
|-----------|---------------|-------|-------------|
| `rate_secretion_immunostimulatory_factor` | 10 | 1/min | Secretion rate of immunostimulatory factor by tumor cells |
| `saturation_density_immunostimulatory_factor` | 1 | - | Saturation density for immunostimulatory factor secretion |

<u>Tumor Cell Glucose Parameters</u>

These parameters only have an effect when `add_glucose` is `true`.

| Parameter | Default Value | Units | Description |
|-----------|---------------|-------|-------------|
| `default_glucose_consumption_tumor_cell` | 0.0007 | 1/min | Baseline glucose consumption rate of tumor cells |
| `glucose_saturation_for_tumor_cell_growth` | 0 | mmol/L | Glucose level above which tumor cells grow at full speed. Below it, growth slows down linearly, which also delays proliferation since cells must grow before dividing |
| `glucose_limit_for_tumor_cell_growth` | 0 | mmol/L | Glucose level below which tumor cells stop growing |
| `glucose_limit_for_death` | 0 | mmol/L | Glucose level below which tumor cells start dying from lack of glucose |
| `glucose_limit_for_death_maximum` | 0 | mmol/L | Glucose level at which the death probability due to lack of glucose is maximum |
| `maximum_death_lack_of_glucose_rate` | 0 | min⁻¹ | Maximum death rate due to lack of glucose |

<u>Tumor Cell Volume Relaxation Rates</u>

| Parameter | Default Value | Units | Description |
|-----------|---------------|-------|-------------|
| `volume_relaxation_rate_alive_tumor_cell_cytoplasm` | 0.00216667 | min⁻¹ | Cytoplasm relaxation rate for alive tumor cells |
| `volume_relaxation_rate_alive_tumor_cell_nucleus` | 0.00366667 | min⁻¹ | Nucleus relaxation rate for alive tumor cells |
| `volume_relaxation_rate_alive_tumor_cell_fluid` | 0.0216667 | min⁻¹ | Fluid relaxation rate for alive tumor cells |
| `volume_relaxation_rate_cytoplasm_necrotic_swelling_tumor_cell` | 5.33333e-05 | min⁻¹ | Cytoplasm relaxation rate during necrotic swelling |
| `volume_relaxation_rate_nucleus_necrotic_swelling_tumor_cell` | 0.000216667 | min⁻¹ | Nucleus relaxation rate during necrotic swelling |
| `volume_relaxation_rate_fluid_necrotic_swelling_tumor_cell` | 0.000833333 | min⁻¹ | Fluid relaxation rate during necrotic swelling |
| `volume_relaxation_rate_cytoplasm_necrotic_lysed_tumor_cell` | 5.33333e-05 | min⁻¹ | Cytoplasm relaxation rate during necrotic lysis |
| `volume_relaxation_rate_nucleus_necrotic_lysed_tumor_cell` | 0.000216667 | min⁻¹ | Nucleus relaxation rate during necrotic lysis |
| `volume_relaxation_rate_fluid_necrotic_lysed_tumor_cell` | 0.000833333 | min⁻¹ | Fluid relaxation rate during necrotic lysis |

<u>Tumor Cell Type Thresholds</u>

| Parameter | Default Value | Units | Description |
|-----------|---------------|-------|-------------|
| `threshold_cancer_cell_type1` | 1.5 | - | Oncoprotein threshold for Type 1 (most aggressive) |
| `threshold_cancer_cell_type2` | 1.0 | - | Oncoprotein threshold for Type 2 |
| `threshold_cancer_cell_type3` | 0.5 | - | Oncoprotein threshold for Type 3 |
| `threshold_cancer_cell_type4` | 0.0 | - | Oncoprotein threshold for Type 4 (least aggressive) |

#### CAR-T Cell Parameters

| Parameter | Default Value | Units | Description |
|-----------|---------------|-------|-------------|
| `average_maximum_time_until_apoptosis_cart` | 12342.86 | minutes | Average CAR-T cell lifespan |
| `default_oxygen_consumption_cart` | 1 | 1/min | Baseline oxygen consumption rate of CAR-T cells |
| `default_glucose_consumption_cart` | 0.0007 | 1/min | Baseline glucose consumption rate of CAR-T cells (only used if `add_glucose` is `true`) |
| `default_volume_new_cart_cell` | 2494 | μm³ | Mean total volume of newly created CAR-T cell |
| `std_volume_new_cart_cell` | 0 | μm³ | Standard deviation of the total volume of newly created CAR-T cells |
| `min_volume_new_cart_cell` | 2494 | μm³ | Minimum total volume of newly created CAR-T cells (sampled volumes are clipped) |
| `max_volume_new_cart_cell` | 2494 | μm³ | Maximum total volume of newly created CAR-T cells (sampled volumes are clipped) |
| `default_fraction_of_volume_for_nucleus_cart_cell` | 0.21652 | - | Nuclear volume fraction of newly created CAR-T cell |
| `default_fraction_fluid_cart_cell` | 0.75 | - | Fraction of the CAR-T cell volume that is fluid |
| `kill_rate_cart` | 0.06667 | min⁻¹ | Rate at which CAR-T attempts to kill attached tumor cell |
| `adhesion_rate_cart` | 0.013 | min⁻¹ | Rate at which CAR-T attempts to attach to tumor cells |
| `max_adhesion_distance_cart` | 18 | μm | Maximum distance for CAR-T to tumor cell attachment |
| `min_adhesion_distance_cart` | 14 | μm | Minimum distance for CAR-T to tumor cell attachment |
| `minimum_distance_from_tumor_to_spawn_cart` | 50 | μm | Minimum distance from tumor boundary to spawn CAR-T cells |

<u>CAR-T Cell Motility Parameters</u>

| Parameter | Default Value | Units | Description |
|-----------|---------------|-------|-------------|
| `persistence_time_cart` | 10 | minutes | Average time before CAR-T changes direction |
| `avg_migration_bias_cart` | 0.5 | - | Mean chemotaxis bias toward immunostimulatory factor (0=random, 1=fully directed) |
| `std_migration_bias_cart` | 0 | - | Standard deviation of the migration bias among CAR-T cells (sampled values are clipped to ≤ 1) |
| `migration_speed_cart` | 5 | μm/min | CAR-T cell migration speed |
| `elastic_constant_cart` | 0.01 | - | Elastic constant for CAR-T cell motility |
| `max_squared_distance_cart_moving_towards_tumor_cell` | (radius_cart + radius_tumor + 1)² | μm² | Squared distance below which a CAR-T cell stops actively moving towards its attached tumor cell (avoids pushing it). Computed from the cell volumes if not given |

### Parameter Configuration

To modify parameters, edit the [`params.json`](params.json) file. Only include parameters you want to change from their default values. Example:

```json
{
  "seed": 1,
  "output_performance_statistics": true,
  "total_minutes_to_simulate": 43200,
  "initial_spherical_tumor_radius": 150.0,
  "treatment": {
    "0": 3957,
    "8": 3957
  }
}
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

./build/CARTopiaX
```

---

## Visualizing Results

Data about tumor growth, different types of cell populations and oxygen and oncoprotein levels are output in `./output/final_data.csv` To visualize plots that give an overview of the results of the simulation you can run the provided python notebook `./CARTopiaX_Simulation_Analysis.ipynb`

To visualize the 3D model of the execution in ParaView use:
```bash
paraview ./output/CARTopiaX/CARTopiaX.pvsm

```

---

## Model Calibration (Bayesian Optimization)

Many parameters of the model cannot be measured directly in the laboratory. [`optimize.py`](abm_calibration/optimize.py) finds the values of selected parameters that make the simulation reproduce some target data (e.g. experimental measurements) as closely as possible.

### How it works

Each simulation is expensive, so trying every combination of parameters (grid search) is not feasible. Instead, the script uses **Bayesian optimization** through [Optuna](https://optuna.org/) and its **Tree-structured Parzen Estimator (TPE)** sampler:

1. Each *trial* proposes a value for every parameter to be calibrated within its range.
2. The simulation is run with those values `NUMBER_MONTE_CARLO` times with different random seeds, and the error against the target data is averaged. This reduces the effect of the stochastic nature of the model.
3. TPE builds a probabilistic model from all previous trials, separating the parameter values that gave good results from the ones that gave bad results. The next trial is sampled where good values are more likely, balancing *exploration* of new regions and *exploitation* of promising ones.
4. After `NUMBER_OF_TRIALS` trials, the parameters with the lowest error are reported.

This way good parameters are usually found with far fewer simulations than with a random or grid search.

### Error modes

The error is selected with the `MODE` variable:

| Mode | Target file | Description |
|------|-------------|-------------|
| `total` | `target_data/final_data.csv` | Mean squared error of a metric of the whole tumor over time (e.g. `average_oxygen_all_cells`, `tumor_radius`) |
| `radius` | `target_data/data_dependent_on_radius_tumor.csv` | Mean squared error of the radial profile of a metric at a fixed minute. Requires `output_information_dependent_on_radius: true` |
| `custom` | - | User-defined error implemented in `compute_custom_error()` |

Target CSVs must have the same format as the files written by the simulation in `output/`, and are compared on the common values of `total_minutes`.

### Configuration

The sections of the script marked with `Change this` are meant to be adapted to each experiment:

- **Settings:** `EXPERIMENT_ID`, `MODE`, `SEED`, `NUMBER_OF_TRIALS`, `NUMBER_MONTE_CARLO`, `TOTAL_MINUTES_TO_SIMULATE` and `RADIUS_PROFILE_MINUTE` (minute of the radial profile compared in `radius` and `custom` modes, which must not be later than `TOTAL_MINUTES_TO_SIMULATE`).
- **`run_ABM()`:** fixed simulation configuration for each run (same keys as `params.json`), together with the parameters being calibrated.
- **`objective()`:** parameters to be calibrated and their search ranges, e.g. `trial.suggest_float("initial_oxygen_level", 30, 40)`.
- **Error functions:** metric to compare.

### Running the calibration

The simulations run inside the Python process through the `cppyy` bundled with BioDynaMo's ROOT, and their output is read in memory instead of through `params.json` and the CSV files. [Build the simulation](#building-the-simulation) first (the script loads `build/libCARTopiaX.so`), then install the dependencies, source BioDynaMo and run the script:
```bash
pip install -r abm_calibration/requirements.txt
source <path-to-biodynamo>/bin/thisbdm.sh
python3 abm_calibration/optimize.py
```

Results are stored in `abm_calibration/experiment_<EXPERIMENT_ID>/`:
- `abm_optuna.db`: SQLite database with all trials. Running again with the same `EXPERIMENT_ID` resumes the study instead of starting from scratch. Setting `NUMBER_OF_TRIALS = 0` just prints the best result stored.
- `optuna_results.csv`: table with the parameters and error of every trial.

During the calibration ParaView export is disabled for each run without modifying `bdm.toml`, and `params.json` is not used. If the process is stopped in the middle of a trial (e.g. by a crash inside a simulation), that trial is marked as failed the next time the study is resumed.

---

## Acknowledgments

This project builds upon the BioDynaMo simulation framework.

> Lukas Breitwieser, Ahmad Hesam, Jean de Montigny, Vasileios Vavourakis, Alexandros Iosif, Jack Jennings, Marcus Kaiser, Marco Manca, Alberto Di Meglio, Zaid Al-Ars, Fons Rademakers, Onur Mutlu, Roman Bauer.
> *BioDynaMo: a modular platform for high-performance agent-based simulation*.
> Bioinformatics, Volume 38, Issue 2, January 2022, Pages 453–460.
> [https://doi.org/10.1093/bioinformatics/btab649](https://doi.org/10.1093/bioinformatics/btab649)

CARTopiaX is based on the mathematical models and solver implementations from the research of
Luciana Melina Luque and collaborators, to replicate its findings:

> Luque, L.M., Carlevaro, C.M., Rodriguez-Lomba, E. et al.
> *In silico study of heterogeneous tumour-derived organoid response to CAR T-cell therapy*.
> Scientific Reports 14, 12307 (2024).
> [https://doi.org/10.1038/s41598-024-63125-5](https://doi.org/10.1038/s41598-024-63125-5)

---

## License

This project is licensed under the Apache License 2.0. See the [LICENSE](LICENSE) file for details.

---
## Author Contact Information
**Author:** Salvador de la Torre Gonzalez  
You can check my profile at Princeton University's [Compiler Research Team](https://compiler-research.org/team/SalvadordelaTorreGonzalez). 
Do not hesitate to reach out in case of having questions, email: *delatorregonzalezsalvador at gmail.com*

**Coauthor:** [Luciana Melina Luque](https://www.lmluque.com/)
