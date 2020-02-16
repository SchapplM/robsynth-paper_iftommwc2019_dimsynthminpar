# Exploiting Dynamics Parameter Linearity for Design Optimization in Combined Structural and Dimensional Robot Synthesis


Authors:  
Moritz Schappler, Svenja Tappe, Prof. Dr.-Ing. Tobias Ortmaier

This code is supplemental material to the following submission to IFToMM World Congress 2019:

```bib
@Conference{SchapplerTapOrt2019b,
  author    = {Schappler, M. and Tappe, S. and Ortmaier, T.},
  title     = {Exploiting Dynamics Parameter Linearity for Design Optimization in Combined Structural and Dimensional Robot Synthesis},
  booktitle = {Proc. of the 15th IFToMM World Congress},
  year      = {2019},
  month     = {October},
  doi       = {10.1007/978-3-030-20131-9_193},
}
```

## Abstract

In the design optimization of robot manipulators regarding drive train and link geometries the dynamics equations have to be evaluated repeatedly.
The method proposed in this paper reduces the computational effort in the dynamics evaluations by using the property of parameter linearity of the dynamics equations.
The combined structural and dimensional synthesis of robot manipulators is adapted in a set of hierarchical optimization loops to exploit this dynamics property.
By this means a reduction of computation time for the inverse dynamics in the synthesis of up to factor three is possible.

## Contents

This repository contains Matlab scripts and Latex code to reproduce all figures of the paper and to check the given number of calculations for the different robot models.  
Matlab scripts:

  * `calc/dimsynth_timing_example.m` creates data for Fig. 4
  * `paper/figures/dyn_plin_fcn_comp.m` creates Fig. 4 from the data
  * `proof` contains matlab files for dynamics with the number of calculations in Tab. 1 of the paper and and Excel sheet assembling the numbers.
  * `calc/regressor_properties.m` shows properties of row E. of Tab. 1 in the paper

## Prerequisites

For the optimization and the figure generation to work, the following steps have to be performed:

1. Install Matlab with (tested with R2019a)
2. Set up the Matlab mex compiler, if not already done
3. Download Matlab dependencies containing general robot models and Matlab functions:
  * https://github.com/SchapplM/robotics-toolbox (tested with rev. bcb2e77; 2020-02-15)
  * https://github.com/SchapplM/robsynth-serroblib (tested with rev. 03793751f; 2020-02-15)
  * https://github.com/SchapplM/robsynth-modelgen (tested with rev. 0f11fd8; 2020-02-15)
  * https://github.com/SchapplM/matlab_toolbox (tested with rev. 366f149; 2020-02-15)
  * https://github.com/SchapplM/robotics-dep-ext (tested with rev. a437871; 2020-02-15)
4. Set up the toolboxes to work and modify the following files according to the instructions
  * robsynth-serroblib/maplerepo_path.m.template (has to link to robsynth-modelgen)
  * robsynth-modelgen/robotics_repo_path.template (has to link to robsynth-robotics_toolbox)
5. Run the path initialization scripts (`..._path_init.m`) of all downloaded software repos in Matlab. Put the run commands in the startup.m file of Matlab for easier use.

