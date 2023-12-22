# Certifiably Correct Range-Aided SLAM (CORA)

This is a (no longer supported) MATLAB implementation build on [manopt](https://www.manopt.org/), primarily as a demo piece for extension and comparison.

We recommend you look at our [official C++ implementation](https://github.com/MarineRoboticsGroup/cora) for a performant and actively maintained library.

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

Any feedback, issues, or contributions are welcome and encouraged! Although the development on this is mostly dead, we'll do the best we can to provide support for hopeful users.

Enjoy!
