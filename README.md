# Certifiably Correct Range-Aided SLAM (CORA)

This is the official repository for the paper ["Certifiably Correct Range-Aided SLAM"](https://arxiv.org/abs/2302.11614) by
[Alan Papalia](https://alanpapalia.github.io), Andrew Fishberg, Brendan O'Neill, [Jonathan P. How](https://www.mit.edu/~jhow/),
[David M. Rosen](https://david-m-rosen.github.io/) and [John J. Leonard](https://meche.mit.edu/people/faculty/JLEONARD@MIT.EDU).

For now we provide a simple MATLAB implementation build on [manopt](https://www.manopt.org/), primarily as a demo piece for extension and comparison.

## Examples

Here is CORA in action on the Plaza1 and Plaza2 experiments from the
[CMU Navigating with Ranging Radios Dataset](https://infoscience.epfl.ch/record/283435). 
In these animations the solver is initialized randomly to
emphasize that CORA's performance is agnostic to the provided initialization.

| Plaza2      | Plaza2     |
| ------------| -----------|
| ![Plaza2](media/pz1_grey_solid_inf_loop.gif) | ![Plaza2](media/pz2_grey_solid_inf_loop.gif) |


## Running Experiments

You can recreate all of the experiments from our paper using the data and scripts in [this repository](https://github.com/MarineRoboticsGroup/cora-experiments/tree/main).

## Feedback, Support, and Contributions

Any feedback, issues, or contributions are welcome and encouraged! We'll do the best we can to provide support for hopeful users.

A C++ implementation is a hopeful future contribution. Please reach out if you are interested in contributing or have use for this, as it may affect how we prioritize development!

Enjoy!
