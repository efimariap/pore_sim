# pore_sim
Pore Simulation and Nearest-Neighbor Analysis

This repository contains MATLAB code for generating periodic pore surfaces, detecting pore features using thresholding, and computing spatial statistics such as the Nearest Neighbor Index (NNI). The workflow simulates the surfaces, extracts object centroids, and evaluates structural organization across multiple parameters and iterations.

---

## Simulation Parameters

Key parameters inside the script include:

* **Nx, Ny** — size of the generated surface
* **meandgpx, meandgpy** — spacing between mound centers
* **meanwgpx, meanwgpy** — mound width
* **meanhgp** — mound height
* **sdhgp_values** — standard deviation of mound heights
* **thresholds** — threshold sweep from 0.05 to 0.95
* **numIterations** — repeats of each condition


This project is released under the **MIT License**.
See the `LICENSE` file for details.


---

