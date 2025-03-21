import subprocess
import os
import numpy as np
import pandas as pd
import starfile

from eulerangles import euler2matrix, invert_rotation_matrices, matrix2euler
from pathlib import Path


def get_rotmat(s):
    """
    Generate an (N,3,3) array of a rotation matrix for a given icosahedral rotation matrix
    """
    symOps = open("symops.txt", "w")
    subprocess.run(["relion_refine", "--print_symmetry_ops", "--sym", s], stdout=symOps)
    a = []
    with open("symops.txt") as f:
        for lines in f:
            lines = lines.strip()
            if "R" in lines:
                continue
            if "+" in lines:
                continue
            if "Euler" in lines:
                continue
            else:
                a.append(lines.split())
    sym_rot_mat = np.array(a, dtype="float").reshape(-1, 3, 3)
    os.remove("symops.txt")
    return sym_rot_mat


def star_downgrade(star_in):
    """
    Downgrade the relion 3.1.1 starfile to relion 3.0.
    """
    cmd = "relion_star_downgrade -s".split()
    cmd.insert(2, star_in)
    print(cmd)
    subprocess.run(cmd)


def anySame(x, y, tol=0.001):
    """
    Filter matrices based on the difference. This became necessary due to
    interpolation errors introduced during the dot multiplication
    """
    for i in x:
        diff = np.abs(np.array(i) - np.array(y))
        if (diff < tol).all():
            return False
    return True


def get_subparticle_vector(sym, subsym):
    """
    Returns a subparticle basis vector based on the selection of a global icosahedral symmetry and
    the sub symmetry
    e.g. t = get_subparticle_vector("I2", "C5")
    """
    sym = sym.upper()
    subsym  = subsym.upper()

    i1five = np.array(
        [0, np.sin(np.radians(31.7175)), np.cos(np.radians(31.7175))]
    ).reshape(3, 1)
    i1three = np.array(
        [np.sin(np.radians(20.91)), 0, np.cos(np.radians(20.91))]
    ).reshape(3, 1)
    i1two = np.array([0, 0, 1]).reshape(3, 1)

    i2five = np.array(
        [np.sin(np.radians(31.7175)), 0, np.cos(np.radians(31.7175))]
    ).reshape(3, 1)
    i2three = np.array(
        [0, np.sin(np.radians(20.91)), np.cos(np.radians(20.91))]
    ).reshape(3, 1)
    i2two = np.array([0, 0, 1]).reshape(3, 1)

    i3five = np.array([0, 0, 1]).reshape(3, 1)
    i3three = np.array(
        [0, np.sin(np.radians(20.91)), np.cos(np.radians(20.91))]
    ).reshape(3, 1)
    i3two = np.array(
        [np.sin(np.radians(31.7175)), 0, np.cos(np.radians(31.7175))]
    ).reshape(3, 1)

    if sym == "I1":
        if subsym == "C5":
            vector = i1five
        elif subsym == "C3":
            vector = i1three
        elif subsym == "C2":
            vector = i1two
    elif sym == "I2":
        if subsym == "C5":
            vector = i2five
        elif subsym == "C3":
            vector = i2three
        elif subsym == "C2":
            vector = i2two
    elif sym == "I3":
        if subsym == "C5":
            vector = i3five
        elif subsym == "C3":
            vector = i3three
        elif subsym == "C2":
            vector = i3two
    return vector

def custom_vector(xyz, box_size):
    '''
    A custom vector can be defined by opening the consensus average in IMOD
    and placing querying 3dmod to report the coords under the cursor. Alternatively,
    a model file can be generated for the point the user is interested in, converted
    to a text file and read in.
    The coords from IMOD are relative to the 0,0,0 from the bottom left of the volume.
    To recenter on the consensus average, read in the coords and subtract half the box size in pixels
    and normalise the vector.
    '''
    half_box = box_size / 2
    point = np.array(xyz.split(','), dtype='float')
    point = point - half_box
    return point / np.linalg.norm(point, ord=2)

def compute_matrix(vector):
    '''
    Lifted from localised reconstuction...
    '''
    if abs(vector[0] < 0.00001) and abs(vector[1] < 0.00001):
        rot = 0.00
        tilt = 0.00
    else:
        rot = np.arctan2(vector[1], vector[0])
        tilt = np.arccos(vector[2])

    psi = 0 #What happens if this changes? Test it out.

    rot1 = float(np.degrees(rot))
    tilt1 = float(np.degrees(tilt))
    psi1 = float(psi)
    eulers = np.array([rot1, tilt1, psi1], dtype=float)
    vector_matrix = euler2matrix(
        eulers, axes="zyz", intrinsic=True, right_handed_rotation=True
    )
    return vector_matrix

def axisAngle2Matrix(axis_angle=np.array([0, 0, 1, 0])):
    """
    https://www.euclideanspace.com/maths/geometry/rotations/conversions/angleToMatrix/index.htm
    Convert angle, axis to rotation matrix
    input is a 1,4 vector where first three elements represent the x,y,z axis
    and the 4th element is the angle in degrees
    """
    angle = np.radians(axis_angle[3])
    c = np.cos(angle)
    s = np.sin(angle)
    t = 1 - c
    x = axis_angle[0]
    y = axis_angle[1]
    z = axis_angle[2]

    turn = np.array(
        [
            [t * x * x + c, t * x * y - z * s, t * x * z + y * s],
            [t * x * y + z * s, t * y * y + c, t * y * z - x * s],
            [t * x * z - y * s, t * y * z + x * s, t * z * z + c],
        ]
    )
    return turn


# get user input
STAR_FILE = Path("Refine3D/I3sym_virion/run_ct21_data.star")
OUTPUT_STAR_FILE = "Refine3D/I3sym_virion/run_ct21_data_harold.star"
ICOS_CONVENTION = "I3"
SUBPARTICLE_VECTOR = get_subparticle_vector(ICOS_CONVENTION, "C5")

# get particle info
star = starfile.read(STAR_FILE)
euler_headings = [f"rlnAngle{axis}" for axis in ("Rot", "Tilt", "Psi")]
coordinate_headings = [f"rlnCoordinate{axis}" for axis in "XY"]
tomogram_heading = "rlnMicrographName"

eulers_particles = star['particles'][euler_headings].to_numpy()
coords_particles = star['particles'][coordinate_headings].to_numpy()
tomograms_particles = star['particles'][tomogram_heading].to_numpy()

# get particle rotation matrices
# these matrices rotate basis vectors in the reference frame of the
# reference to aign with the reference frame of the particle
particle_rotation_matrices = invert_rotation_matrices(
    euler2matrix(
        eulers_particles, axes="zyz", intrinsic=True, right_handed_rotation=True
    )
)

# get icosahedral matrices (60)
icosahedral_matrices = get_rotmat(ICOS_CONVENTION)

# find set of coincedent rotated z-vectors (subparticles)
rotated_z_vectors = (icosahedral_matrices @ SUBPARTICLE_VECTOR).squeeze()

# find subparticle rotation matrices and shifts
already_seen = []
subparticle_rotation_matrices = []

for idx, rm in enumerate(icosahedral_matrices):
    rotated_z_vector = tuple(rotated_z_vectors[idx])
    fresh = anySame(already_seen, rotated_z_vector)
    if fresh:
        already_seen.append(rotated_z_vector)
        subparticle_rotation_matrices.append(rm)

subparticle_rotation_matrices = np.array(subparticle_rotation_matrices)  # (12, 3, 3)
subparticle_orientations = subparticle_rotation_matrices.reshape(
    -1, 1, 3, 3
)  # (m, 1, 3, 3)



# write out a star file with subparticle info
newStar=star.copy()

# this is the dataframe corresponding to particle information as opposed to optics groups etc.
parts=pd.DataFrame(star['particles'])

#create an empty datafram to fill with symmetry expanded particles
newParts=pd.DataFrame(star['particles']).iloc[0:1]

N_sym=np.shape(subparticle_orientations)[0]
N_parts=np.shape(particle_rotation_matrices)[0]

for i in np.arange(N_parts):

    print("particle " +str(i) +" of " +str(N_parts) + " finished. ("+str(100*float(i)/float(N_parts))+"%)" , end = "\r")
    #create a dataframe to hold the current particles expanded subparticles
    expandedParts=pd.DataFrame([parts.iloc[i]]*N_sym)

    #get the orientation of this particle
    thisParticleOri = particle_rotation_matrices[i]
    align_subparticle_on_z = compute_matrix(SUBPARTICLE_VECTOR)
    expandedPartOris = thisParticleOri @ subparticle_orientations @ align_subparticle_on_z
    expandedPartEuls = matrix2euler(invert_rotation_matrices(expandedPartOris.reshape((-1, 3, 3))),
				     axes="zyz",
    				     intrinsic=True,
    				     right_handed_rotation=True)
    expandedParts['rlnAngleRot'] = expandedPartEuls[:,0]
    expandedParts['rlnAngleTilt'] = expandedPartEuls[:,1]
    expandedParts['rlnAnglePsi'] = expandedPartEuls[:,2]

    newParts=newParts.append(expandedParts)

#newParts.reset_index(inplace=True)
newStar['particles']=newParts.iloc[1:]


starfile.write(newStar, OUTPUT_STAR_FILE, overwrite=True)
print(f"Wrote out {N_parts*N_sym} subparticles in {OUTPUT_STAR_FILE}")
