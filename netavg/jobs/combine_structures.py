from os.path import join
from glob import glob
from prody import parsePDB, writePDB


MAX_FRAMES = 1000


def combine_structures(directory_with_pdbs, output_filename):
    search_string = join(directory_with_pdbs, '*.pdb')
    path_list = glob(search_string)

    if len(path_list) > MAX_FRAMES:
        raise RuntimeError('Got %d frames, but only up to %d frames are allowed.')
    else:
        pass

    atom_group = parsePDB(path_list[0])
    for i, path in enumerate(path_list):
        if i == 0:
            continue
        else:
            p = parsePDB(path)
            atom_group.addCoordset(p)
    writePDB(output_filename, atom_group)
