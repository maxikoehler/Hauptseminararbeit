# smib-simulation
A dynamic simulation of a single machine infinite busbar system. Please see the corresponding [Blogpost on Medium](https://medium.com/@georg.kordowich/watts-up-with-dynamic-power-system-simulations-c0f16fc99769)
for more details.

The repository contains two python files to run the simulation of a single machine infinite busbar model.

- `simulation.py` contains the code for the simulation using Euler's method.
- `simulation_even_better.py` contains the code for the simulation using Heun's method.
- `simulation.ipynb` is a jupyter notebook which contains the code for the simulation using Euler's method and is a 
duplicate of `simulation.py` so you can run it directly on binder: [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/georgkordowich/smib-simulation/HEAD?labpath=simulation.ipynb)
