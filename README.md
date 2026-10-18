# Stabilization of Nonlinear Systems with State-Dependent Representation

This repository contains the MATLAB and experimental code accompanying the
paper _Stabilization of Nonlinear Systems with State-Dependent
Representation: From Model-Based to Direct Data-Driven Control_.

The code studies controller synthesis for discrete-time nonlinear systems
written in the state-dependent form

$$
x^{+}=A(x)x+B(x)u.
$$

The examples cover model-based and direct data-driven stabilization,
disturbance robustness, input saturation, region-of-attraction estimation,
quadrotor simulation, and physical experiments.

## Resources

- [Demonstrations and videos](https://nus-core.github.io/assets/standalone/AxxBxu/index.html)
- [Preprint on arXiv](https://arxiv.org/abs/2510.16451)

## Requirements

Depending on the selected script, the MATLAB workflows use:

- MATLAB;
- [CVX](https://cvxr.com/cvx/) for the CVX-based semidefinite programs;
- [YALMIP](https://yalmip.github.io/) for scripts built with
  `sdpvar`, `optimize`, or `solvesos`;
- [MOSEK](https://www.mosek.com/) as the configured SDP/SOS solver.

Several scripts explicitly select MOSEK with `cvx_solver mosek` or through
YALMIP solver settings. To run those scripts unchanged, install MOSEK,
activate a valid license, and make sure it is visible to the corresponding
modeling toolbox. A different compatible solver may be selected where the
formulation and script settings permit it.

The Gazebo and physical-flight components in `SimE/` additionally use
ROS 1, Gazebo, MAVROS, Eigen, and the
[XTDrone environment](https://www.yuque.com/xtdrone).

The associated XTDrone paper is available [this site](https://arxiv.org/abs/2003.09700).

## Quick start

Clone or download the repository, configure CVX or YALMIP as required by
the selected script, and verify the MOSEK installation. Run each script
from its own directory because the data paths are relative to that
directory. For example:

~~~matlab
cd SimB
Model_stabilization_04
~~~

To run the quadrotor simulation example:

~~~matlab
cd SimC
Quadrotor08
~~~

To synthesize the controller from the processed physical data:

~~~matlab
cd SimE
controller_synthesize
~~~

The scripts place their main results in the MATLAB workspace and generate
the corresponding figures. Some examples also read or write files in their
local `Mdata/` or `output/` directories.

## Paper-to-code map

The section and figure references below follow the updated manuscript
[`paper_2510.16451v2.pdf`](SimB/previous/paper_2510.16451v2.pdf).

| Paper location | Main script or directory | Description |
| --- | --- | --- |
| Section IV-A, Table I | `SimB/Data_compare_Axx_02.m` and `SimB/Data_compare_SOS_02.m` | Computational scalability comparison with the SOS baseline |
| Section IV-A, Fig. 1 | `SimB/Data_cancellation_04.m` | ROA estimation comparison |
| Section IV-B | `SimE/` | Physical experiment and supporting Gazebo workflow |
| Appendix B-A, Fig. 3 (left) | `SimB/Model_stabilization_04.m` | Model-based closed-loop phase portrait |
| Appendix B-A, Fig. 3 (right) | `SimB/Model_stabilization_04_r_and_r0.m` | Model-based ROA radius versus prescribed state-region radius |
| Appendix B-A, Fig. 4 | `SimB/Model_stabilization_04.m` | Model-based state response under disturbance |
| Appendix B-A, Fig. 5 | `SimB/Model_saturation_04.m` | Model-based ROA estimation with input saturation |
| Appendix B-B, Fig. 6 (left) | `SimB/Data_stabilization_04.m` | Data-driven closed-loop phase portrait |
| Appendix B-B, Fig. 6 (right) | `SimB/Data_stabilization_04_r_and_r0.m` | Data-driven ROA radius versus prescribed state-region radius |
| Appendix B-B, Fig. 7 | `SimB/Data_stabilization_04.m` | Data-driven state response under disturbance |
| Appendix B-B, Fig. 8 | `SimB/Data_saturation_04.m` | Data-driven ROA estimation with input saturation |

`SimC/Quadrotor08.m` remains available as an additional quadrotor numerical
example; it is not part of the paper-to-code map above.

## Repository structure

| Directory | Purpose |
| --- | --- |
| `SimB/` | Model-based and data-driven simulations and comparative studies |
| `SimC/` | Quadrotor simulation |
| `SimE/` | Gazebo and physical-experiment workflow |

## Gazebo and physical-experiment data

The principal `SimE/` files are:

- `data_physical_exp_collect_04.m`: imports and aligns the physical CSV
  measurements;
- `controller_synthesize.m`: constructs the data matrices and synthesizes
  the feedback controller with CVX;
- `ros.AxxBxu.cpp`: reference ROS/Gazebo node source.

The repository does not contain a complete catkin package for
`ros.AxxBxu.cpp`. The source includes project-specific `my_control` headers and
must be integrated into an existing ROS/XTDrone catkin workspace with its
ROS message and build dependencies configured.

Processed CSV and MAT data used by the supplied scripts are included under
`SimE/Gazebo_data_csv/`, `SimE/physical_data_csv/`,
`SimE/Gazebo_data_mat/`, and `SimE/physical_data_mat/`. The original
ROS bag files are not tracked because of their size and can be downloaded
from the
[external data folder](https://nusu-my.sharepoint.com/:f:/g/personal/li_ld_nus_edu_sg/Esiid6VxqM9AnskADErBVK8BECar4wHoJvvaXhQw8TRjeQ?e=KSHiYJ).
By default, `controller_synthesize.m` loads the processed physical data
`SimE/physical_data_mat/data_exp_09A.mat`; the Gazebo MAT-file alternative
is retained as a commented option in the script.

## Reproducibility notes

- Some examples generate inputs, initial conditions, noise, or disturbances
  at run time. Those results may vary unless the script fixes the
  random-number seed; examples that load supplied MAT files reuse fixed data.
- Run scripts from the directory indicated above so that relative data
  paths resolve correctly.
- Solver status and LMI feasibility should be checked before interpreting
  a returned controller or region-of-attraction estimate.
- Generated outputs may overwrite existing files in an example's
  `output/` directory; preserve previous results separately when needed.

## Citation

If this repository is useful in your work, please cite:

> L. Li, R. Huang, and L. Zhao, “Stabilization of Nonlinear Systems with
> State-Dependent Representation: From Model-Based to Direct Data-Driven
> Control,” _arXiv preprint arXiv:2510.16451_, 2025.
