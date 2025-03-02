# FLIP-Fluid Simulation (C++ Rewrite)

This project is a C++ rewrite of the [FLIP-Fluid simulation](https://github.com/matthias-research/pages/blob/master/tenMinutePhysics/18-flip.html) from Matthias Müller's _Ten Minute Physics_ series. The FLIP (Fluid-Implicit Particle) method is a hybrid approach that combines particle-based and grid-based techniques to simulate realistic fluid dynamics efficiently.

## Features

- Fully implemented in C++ for performance and flexibility
- Hybrid particle-grid solver for smooth fluid behavior
- Supports various boundary conditions
- Optimized for real-time simulations

## Installation & Usage

1. Clone the repository:
   ```bash
    git clone https://github.com/burakssen/flip-fluid-simulation
    cd flip-fluid-simulation
   ```
2. Build using CMake
   ```bash
    mkdir build && cd build
    cmake ..
    make
   ```
3. Run the simulation

   ```bash
    ./flip
   ```

# References

- Original implementation: [Matthias Müller's FLIP simulation](https://github.com/matthias-research/pages/blob/master/tenMinutePhysics/18-flip.html)
- FLIP Fluids: Bridging Particles and Grids
