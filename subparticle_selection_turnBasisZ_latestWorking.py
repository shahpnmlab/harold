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

    if sym == "I1" or sym == "i1":
        if subsym == "C5" or subsym == "c5":
            vector = i1five
        elif subsym == "C3" or subsym == "c3":
            vector = i1three
        elif subsym == "C2" or subsym == "c2":
            vector = i1two
    elif sym == "I2" or sym == "i2":
        if subsym == "C5" or subsym == "c5":
            vector = i2five
        elif subsym == "C3" or subsym == "c3":
            vector = i2three
        elif subsym == "C2" or subsym == "c2":
            vector = i2two
    elif sym == "I3" or sym == "i3":
        if subsym == "C5" or subsym == "c5":
            vector = i3five
        elif subsym == "C3" or subsym == "c3":
            vector = i3three
        elif subsym == "C2" or subsym == "c2":
            vector = i3two
    return vector


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


def compute_matrix(vector):
    if abs(vector[0] < 0.00001) and abs(vector[1] < 0.00001):
        rot = 0.00
        tilt = 0.00
    else:
        rot = np.arctan2(vector[1], vector[0])
        tilt = np.arccos(vector[2])

    psi = 0

    rot1 = np.degrees(rot)
    tilt1 = np.degrees(tilt)
    psi1 = np.degrees(psi)
    eulers = np.array([rot1, tilt1, psi1], dtype=float)
    vector_matrix = euler2matrix(
        eulers, axes="zyz", intrinsic=True, right_handed_rotation=True
    )
    return vector_matrix


# get user input
#STARFILE311 = "run_data.star"
#star_downgrade(STARFILE311)
STAR_FILE = Path("run_data.star")
OUTPUT_STAR_FILE = f"{STAR_FILE.parent / STAR_FILE.stem}_subparticles.star"
ICOS_CONVENTION = "I2"
ICOSAHEDRAL_RADIUS_PX = 32  # Pixels
SUBPARTICLE_VECTOR = get_subparticle_vector(ICOS_CONVENTION, "C5")
AXISANGLE = np.array([0, 1, 0, 31.7175])

# get particle info
star = starfile.read(STAR_FILE)
euler_headings = [f"rlnAngle{axis}" for axis in ("Rot", "Tilt", "Psi")]
coordinate_headings = [f"rlnCoordinate{axis}" for axis in "XYZ"]
tomogram_heading = "rlnMicrographName"

eulers_particles = star[euler_headings].to_numpy()
coords_particles = star[coordinate_headings].to_numpy()
tomograms_particles = star[tomogram_heading].to_numpy()

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
icosahedral_matrices = invert_rotation_matrices(icosahedral_matrices)

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

subparticle_z_vectors = subparticle_rotation_matrices @ SUBPARTICLE_VECTOR  # (12, 3, 1)

subparticle_shifts = ICOSAHEDRAL_RADIUS_PX * subparticle_z_vectors  # (12, 3)
subparticle_shifts = subparticle_shifts.reshape(-1, 1, 3, 1)  # (m, 1, 3, 1)

### apply transformations onto particles
particle_positions = coords_particles.reshape((-1, 3, 1))  # (n, 3, 1)
particle_orientations = particle_rotation_matrices  # (n, 3, 3)

# orient the shift vectors according to particle orientation
oriented_shifts = particle_orientations @ subparticle_shifts  # (m, n, 3, 1)

# apply the transformed shifts onto particle positions
subparticle_positions = particle_positions + oriented_shifts  # (m, n, 3, 1)

# calculate the new orientations by composing rotations
# (n, 3, 3) @ (m, 1, 3, 3) -> (m, n, 3, 3)

align_subparticle_on_z = axisAngle2Matrix(AXISANGLE)

orientations = particle_orientations @ subparticle_orientations @ align_subparticle_on_z


output_positions = subparticle_positions.reshape((-1, 3))
output_orientations = orientations.reshape((-1, 3, 3))

output_oriented_z = output_orientations @ SUBPARTICLE_VECTOR

# get eulers from rotation matrices for subparticles

output_eulers = matrix2euler(
    invert_rotation_matrices(output_orientations),
    axes="zyz",
    intrinsic=True,
    right_handed_rotation=True,
)

# output_eulers[:,[0,2]] = -output_eulers[:,[0,2]] # Negate the rot and psi angles

# generate (m, n) set of tomogram names
output_tomograms = np.empty(subparticle_positions.shape[:2], dtype=object)
output_tomograms[:, :] = tomograms_particles

# reshape output tomograms to (m * n)
output_tomograms = output_tomograms.reshape(-1, 1)

# write out a star file with subparticle info
df = pd.DataFrame(data=output_positions, columns=coordinate_headings)
df[tomogram_heading] = output_tomograms
df[euler_headings] = output_eulers
df = df.sort_values(by=['rlnMicrographName'])

starfile.write(df, OUTPUT_STAR_FILE, overwrite=True)
print(f"Wrote out {output_tomograms.shape[0]} subparticles in {OUTPUT_STAR_FILE}")

# Visualise for debug
# import napari
# napari_vectors = np.zeros((output_positions.shape[0], 2, 3))
# napari_vectors[:, 0, :] = output_positions
# napari_vectors[:, 1, :] = output_oriented_z.squeeze()

# v = napari.Viewer()
# v.add_points(output_positions)
# v.add_vectors(napari_vectors)
