# Memristive Adaptive Neuron Spatial Optimisation Automata

![Language](https://img.shields.io/badge/Language-MATLAB-orange.svg)
![License](https://img.shields.io/badge/License-GPL%202.0-blue.svg)

## 📌 Overview

**Memristive_Adaptive_Neuron_Spatial_Optimisation_Automata** is a MATLAB-based repository focused on simulating a bio-mimetic Neuron Cell using $SiO_2$-based adaptive Conductive Bridging Random Access Memories (CBRAMs). The project further extends this foundational hardware concept into an emergent optimal transport route network through a continuous Cellular Automata (CA) application.

This repository explores the intersection of neuromorphic engineering, memristive device modeling, and spatial optimization.

## 🗂️ Project Structure

The repository is organized into four primary modules:

* **`/CBRAM Model`**
  Contains the core mathematical and behavioral models for the $SiO_2$-based adaptive CBRAMs. These models simulate the ion migration and filament formation responsible for the device's memristive properties.
* **`/Memristive_Adaptive_Neuron_Spatial_Optimisation_Local`**
  Implements the continuous Cellular Automata (CA) framework. This module is responsible for the spatial optimization and routing algorithms, demonstrating how interconnected memristive "neurons" can dynamically find optimal transport routes.
* **`/Non-Linear Analysis`**
  Provides analytical scripts and tools to evaluate the non-linear dynamics of the bio-mimetic neuron cells and the resulting emergent network behaviors. 
* **`/Simulink`**
  Contains Simulink models (`.slx`) for block-level simulation and visual representation of the memristive neuron circuits and their systemic interactions.

## 🚀 Prerequisites

To run the models and simulations in this repository, you will need:
* **MATLAB** (R2021a or newer recommended)
* **Simulink** (for running files in the `/Simulink` directory)
* Relevant MATLAB toolboxes (e.g., Control System Toolbox, Simscape—depending on specific script requirements).

## 🛠️ Usage

1. **Clone the repository:**
   ```bash
   git clone [https://github.com/chipher-zero/Memristive_Adaptive_Neuron_Spatial_Optimisation_Automata.git](https://github.com/chipher-zero/Memristive_Adaptive_Neuron_Spatial_Optimisation_Automata.git)
   ```
2. **Open in MATLAB:**
   Navigate to the cloned directory within the MATLAB environment.
3. **Run Models:**
   * To explore the memristive properties, run the initialization scripts within the `/CBRAM Model` directory.
   * To view the system-level simulation, open the corresponding `.slx` files in the `/Simulink` folder.
   * To observe the spatial optimization cellular automata, execute the main run scripts located in the `/Memristive_Adaptive_Neuron_Spatial_Optimisation_Local` directory.

## 📄 License

This project is licensed under the **GPL-2.0 License**. See the [LICENSE](./LICENSE) file for more details.