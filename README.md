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

To use Harold, open the script and edit the following parameters
```
STAR_FILE = Path("Refine3D/I3sym_virion/run_ct21_data.star")
OUTPUT_STAR_FILE = "Refine3D/I3sym_virion/run_ct21_data_harold.star"
ICOS_CONVENTION = "I3"
SUBPARTICLE_VECTOR = get_subparticle_vector(ICOS_CONVENTION, "C5")
```

Note that while Harold is agnostic to the underlying symmetries (as it explicitly depends on Relion functions), there are useful shortcuts only for dealing with I, O, and T symmetries. To obtain subparticles with redundancy filtration for D- and C-symmetries, custom vectors need to be provided.
