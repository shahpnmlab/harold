# harold
USAGE:
harold.py <particles.star> <global symmetry> <subsymmetry>  <distance to subparticle>
For example if you want to get all the pentons of particles refined in I2 convention:
harold.py run_data.star I2 C5 128
harold.py <particles.star> <global symmetry>  <custom vector with box size> <distance to subparticle>
On the other hand if the sub-particle is not coincident with one of
the icosahedral symmetries, one can provide a string in which the
first three numbers are the XYZ coords of the subparticle as determined
in a visualisation program such as 3dmod and the fourth number is
the size of the box enclsoing the consensus volume.
For example: harold.py run_data.star I2 "100,86,51,178" 128
