# Harold
Harold is a small python implementation for symmetry expansion and redundancy filteration.

## Requirements
To run harold make sure you have installed the following dependencies via pip
- numpy 
- pandas 
- starfile
- eulerangles

preferably its own conda environment.
You also need to have `relion` binaries defined in your `$PATH` environment variable.

## Guide
### SPA (2D particles)
To use harold the basic syntax is as follows :- 
`harold2d.py <particles.star> <global symmetry> <subsymmetry>`
For example if you want to get all the pentons of particles refined in the icosahedral I2 convention, you can issue the following command -
`harold2d.py run_data.star I2 C5`
Note that while harold is agnotic to the underlying symmetries (as it explictly depends on relion functions) there are usefule shortcuts only for dealing with I, O and T symmetries. To obtain subparticles with redundancy filteration for D- and C-symmetries, custom vectors need to be provided and this is not currently implemented for 2D data.

### STA (3D particles / subtomo data)
You can also provide a custom vector to harold2d to define a subparticle arbitarily in the symmetry convetion of your reconstruction.
`harold.py <particles.star> <global symmetry>  <custom vector with box size>` 
If you provide a a custom vector which does not exactly coincide with a subsymmetry axis then you will end up with a fully exapnded set of orientations. The custom vector can be determined in a 3d map visualisation programs such as `3dmod` a package within the the IMOD software suite (highly recommended). The custom vector is provided as a 4 element comma separated string. For example -
`harold.py run_data.star I2 "100,86,51"`
