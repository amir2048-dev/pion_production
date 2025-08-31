pion_production

A Geant4-based simulation framework for studying pion production and transport in matter.
Developed as part of a Master’s research project, this code models interactions of high-energy particles with a converter, optional magnetic fields, and detector geometries.

Features

Modular Geant4 simulation with separate src/ and include/ directories.

Supports single-thread (ST) and multi-thread (MT) builds.

Configurable geometry: world, absorber, and detector volumes.

Optional uniform magnetic field with user-defined strength and integration settings.

Example macros for running simple test simulations and visualization.

Build Instructions
Prerequisites

Linux (tested on Ubuntu inside VirtualBox VM)

Geant4
 (compiled with data libraries and visualization)

CMake ≥ 3.19

g++ or clang++

Build (Single-threaded)
source /path/to/geant4-install/bin/geant4.sh   # or your setup script
cmake -S . -B build -DGEANT4_USE_MULTITHREADED=OFF
cmake --build build -j8

Build (Multi-threaded)
source /path/to/geant4-install_mt/bin/geant4.sh
cmake -S . -B buildmt -DGEANT4_USE_MULTITHREADED=ON
cmake --build buildmt -j8

Running the Simulation

Run from the build directory with a macro:

./build/sim macro/run1.mac


Or for the MT build:

./buildmt/sim macro/run1.mac


If no macro is provided, the Qt visual interface will launch.

Macros

run1.mac – simple run configuration (single particle, minimal geometry).

vis.mac – visualization setup.

You can edit macros under the macro/ directory, then run them directly. They are also copied into the build directory automatically.

Notes on Fields & MT

The magnetic field is defined in ConstructSDandField() to ensure MT safety.

Do not attach a global FieldManager to the world volume in Construct() — Geant4 manages this automatically.

To change geometry or field values after initialization, run:

/run/beamOn 0
/run/reinitializeGeometry
/run/initialize

Git Workflow (for contributors / future you)

Commit locally as often as you like in VS Code.

Push (git push) to GitHub regularly for backup.

Branch for new features (e.g., feature/magnetic-field).

.gitignore excludes build/, buildmt/, data, and VS Code artifacts.

Status

This project is actively developed for academic research. Current focus:

Validating pion production in converter materials.

Extending to include magnetic field transport.

Adding detectors and noise modeling in later stages.

Amir Eyal, august-30, 2025

