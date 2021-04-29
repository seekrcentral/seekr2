"""

"""

import math
import itertools
import collections
import sys

import mdtraj
import numpy as np

def get_box_vector_dimensions(trajectory):
    """ """
    box_vectors = trajectory.unitcell_vectors[0]
    a, b, c = mdtraj.utils.box_vectors_to_lengths_and_angles(
        *box_vectors
    )[:3]
    return np.array([a, b, c])

def apply_pbc_to_coordinates(
        box_vector_dimensions, atom_coordinates, box_vectors,
        coordinate_components):
    """ """
    if coordinate_components:
        box_vector_component = coordinate_components[-1]
        for coordinate_component in coordinate_components:
            scale_factor = math.floor(
                atom_coordinates[coordinate_component] \
                / box_vector_dimensions[coordinate_component]
        )
            atom_coordinates[coordinate_component] -= \
                scale_factor * box_vectors[
                    box_vector_component
                ][
                    coordinate_component
                ]    
        coordinate_components.pop()
        apply_pbc_to_coordinates(
            box_vector_dimensions,
            atom_coordinates,
            box_vectors,
            coordinate_components
        )
    return atom_coordinates

def apply_pbc(trajectory):
    """
    
    """
    box_vector_dimensions = get_box_vector_dimensions(trajectory)
    box_vectors = trajectory.unitcell_vectors[0]
    coordinates = get_coordinates(trajectory)
    for atom_coordinates in coordinates:
        coordinate_components = collections.deque([0, 1, 2])
        apply_pbc_to_coordinates(
            box_vector_dimensions, 
            atom_coordinates, 
            box_vectors,
            coordinate_components
        )
    return trajectory

def get_coordinates(trajectory):
    """ """
    return [atom_coordinates for atom_coordinates in trajectory.xyz[0]]

def coordinate_bounds(coordinates):
    """ """
    bounds = collections.namedtuple("bounds", ["xmin", "xmax", "ymin", "ymax", 
                                   "zmin", "zmax"])
    return bounds(*[bound for upper_lower in [(
        lambda coordinate: [min(coordinate), max(coordinate)]) \
        (coordinate) for coordinate in list(zip(*coordinates))] \
        for bound in upper_lower]) 
    return

def construct_partition(minimum, maximum, num_subintervals):
    """
    
    """
    delta = (abs(minimum) + abs(maximum)) / num_subintervals
    subinterval_bounds =  np.array([
        minimum + index * delta for index
        in range(num_subintervals + 1)
    ])
    return np.array([
        np.array([
            subinterval_bounds[index], 
            subinterval_bounds[index + 1]
        ]) 
        for index in range(num_subintervals)
    ])
    return 

def construct_cubic_partition(coordinate_bounds, num_cubes):
    """
    
    """
    num_component_subintervals = math.ceil(pow(num_cubes, 1 / 3))
    xpart = construct_partition(
        coordinate_bounds.xmin, 
        coordinate_bounds.xmax, 
        num_component_subintervals
    )
    ypart = construct_partition(
        coordinate_bounds.ymin, 
        coordinate_bounds.ymax, 
        num_component_subintervals
    )
    zpart = construct_partition(
        coordinate_bounds.zmin, 
        coordinate_bounds.zmax, 
        num_component_subintervals
    )
    return np.array(list(itertools.product(
        *np.array([xpart, ypart, zpart])
    )))

def in_partition(val, partition_bounds):
    """ """
    lower_bound, upper_bound = partition_bounds 
    if val >= lower_bound:
        if val <= upper_bound:
            return True
    return False

def atom_in_cube(atom_coordinates, cube):
    """ """
    x, y, z = atom_coordinates
    x_bounds, y_bounds, z_bounds = cube
    in_cube = sum(
        [
            in_partition(x, x_bounds),
            in_partition(y, y_bounds),
            in_partition(z, z_bounds)
        ]
    )
    if in_cube == 3:
        return True
    return False

def is_water(atom):
    """ """
    if atom.residue.name[:3] in ["HOH", "WAT", "H2O"]:
        return True
    return False

def num_atoms_per_cube(
        cubic_partition, coordinates, atoms, non_water_cube_indices=True):
    """ """
    non_water_cube_indices_set = set()
    num_atoms_per_cube_dict = {
        cube_index : 0 for cube_index in range(len(cubic_partition))
    }
    for index, atom_coordinates in enumerate(coordinates):
        for cube_index, cube in enumerate(cubic_partition):
            if atom_in_cube(atom_coordinates, cube):
                if non_water_cube_indices:
                    if not is_water(atoms[index]):
                        non_water_cube_indices_set.add(cube_index)    
                num_atoms_per_cube_dict[cube_index] += 1    
    return num_atoms_per_cube_dict, non_water_cube_indices_set

def atom_density_per_cube(
        cubic_partition, coordinates, atoms, water_only_cube_densities=False):
    """ """
    atoms_per_cube, non_water_cube_indicies = num_atoms_per_cube(
        cubic_partition, coordinates,
        atoms, water_only_cube_densities
    )
    volumes = np.array([
        np.prod([ 
            math.fabs(cube[0][0] - cube[0][1]), 
            math.fabs(cube[1][0] - cube[1][1]), 
            math.fabs(cube[2][0] - cube[2][1])
        ])
        for cube in cubic_partition
    ])
    print("volumes:", volumes)
    atom_density_per_cube_dict = {
        cube_index : num_atoms / volumes[cube_index]
        for cube_index, num_atoms in atoms_per_cube.items()
    }
    print("atom_density_per_cube_dict:", atom_density_per_cube_dict)
    if water_only_cube_densities:
        water_only_cube_densities_dict = {
            cube_index : num_atoms / volumes[cube_index]
            for cube_index, num_atoms in atoms_per_cube.items()
                if cube_index not in non_water_cube_indicies
        }
        return (
            atom_density_per_cube_dict,
            water_only_cube_densities_dict
        )
    return atom_density_per_cube_dict

def mean_atom_density_per_cube(atom_density_per_cube_dict):
    """ """
    return np.mean(np.array(list(atom_density_per_cube_dict.values())))

def check_water_box_for_bubbles(trajectory, num_cubes):
    """
    
    """
    trajectory.center_coordinates
    topology = trajectory.topology
    apply_pbc(trajectory)
    coordinates = get_coordinates(trajectory)
    bounds = coordinate_bounds(coordinates)
    cubic_partition = construct_cubic_partition(bounds, num_cubes)
    atoms = list(topology.atoms)
    atom_density, water_only_density = atom_density_per_cube(
        cubic_partition, coordinates, atoms)
    mean_atom_density = mean_atom_density_per_cube(atom_density) 
    mean_water_only_density = mean_atom_density_per_cube(
        water_only_density
    ) 
    delta_mean = abs(mean_atom_density - mean_water_only_density)
    if delta_mean >= 0.1* mean_atom_density: 
        box_has_a_bubble = True
    else:
        box_has_a_bubble = False 
    return box_has_a_bubble #, mean_atom_density, mean_water_only_density

if __name__ == "__main__":
    num_cubes = 64
    pdb_file = sys.argv[1]
    trajectory = mdtraj.load(pdb_file)
    has_bubble = check_water_box_for_bubbles(trajectory, num_cubes)
    if has_bubble:
        print("Bubble detected in structure.")
    else:
        print("No bubbles detected in structure.")