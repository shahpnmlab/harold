# Harold

Harold is a Python implementation for symmetry expansion and redundancy filtration.

## Requirements

To run Harold, make sure you have installed the following dependencies via pip:

-   numpy
-   pandas
-   starfile
-   eulerangles

It is recommended to use its own conda environment. You also need to have the `relion` binaries defined in your `$PATH` environment variable. To use it either copy the script file into your project directory or place it in a location that is defined in your `$PATH` environment variable.

## Guide

### SPA (2D particles)

To use Harold, the basic syntax is as follows:

`python harold2d.py <particles.star> <global symmetry> <subsymmetry>`

For example, if you want to get all the pentons of particles refined in the icosahedral I2 convention, you can issue the following command:

`python harold2d.py run_data.star I2 C5`

Note that while Harold is agnostic to the underlying symmetries (as it explicitly depends on Relion functions), there are useful shortcuts only for dealing with I, O, and T symmetries. To obtain subparticles with redundancy filtration for D- and C-symmetries, custom vectors need to be provided, and this is not currently implemented for 2D data.
