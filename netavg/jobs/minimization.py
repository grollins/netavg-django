#!/usr/bin/env python
# encoding: utf-8


import argparse
import prody
import os
import shutil
import subprocess
import numpy
from os.path import join


GMX_PATH = '/usr/local/gromacs/bin/'


mdp_string = '''
define = -DPOSRES
integrator = {integrator}
nsteps = 1000
emtol = 1
nstlist = 1
coulombtype = Cut-off
vdwtype = Cut-off
ns_type = simple
rlist = 1.8
rcoulomb = 1.8
rvdw = 1.8
pbc = xyz
implicit_solvent = GBSA
gb_algorithm = OBC
sa_algorithm = ACE-approximation
rgbradii = 1.8
;nstxout = 1
'''

def parse_args():
    parser = argparse.ArgumentParser(description='Generate trajectory with gaussian flucutations.')
    parser.add_argument('pdb_file', metavar='INPUT_PDB_FILE', help='path to input pdb file')
    parser.add_argument('trajectory', metavar='TRAJECTORY', help='path to input trajectory')
    parser.add_argument('out_file', metavar='OUTPUT_PDB_FILE', help='path to input pdb file')
    parser.add_argument('--pos_res_k', type=float, default=1000.)
    args = parser.parse_args()

    return (args.pdb_file, args.trajectory, args.out_file, args.pos_res_k)


class Minimizer(object):
    def __init__(self, input_pdb_filename, trajectory_filename, force_const):
        self.input_pdb = self._load_pdb(input_pdb_filename)
        self.trajectory = self._load_pdb(trajectory_filename)

    def _load_pdb(self, in_file):
        protein = prody.parsePDB(in_file)
        return protein

    def _get_closest_frame(self):
        output = prody.AtomGroup('Cartesian average coordinates')
        output.setCoords( self.trajectory.getCoords() )
        output.setNames( self.trajectory.getNames() )
        output.setResnums( self.trajectory.getResnums() )
        output.setResnames( self.trajectory.getResnames() )

        ensemble = prody.PDBEnsemble(self.trajectory)
        ensemble.setCoords(self.input_pdb)
        ensemble.superpose()
        rmsds = ensemble.getRMSDs()
        min_index = numpy.argmin(rmsds)

        output.setCoords( ensemble.getCoordsets(min_index) )
        return output

    def _create_no_h_file(self, output_stream):
        # make the index file
        cmd = join(GMX_PATH, 'make_ndx')
        cmd += ' -f min_round_2.gro -o no_h.ndx'
        p1 = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE,
                              stdout=output_stream, stderr=output_stream)
        p1.communicate('q\n')

        # run editconf
        edit_cmd = join(GMX_PATH, 'editconf')
        edit_cmd += ' -f min_round_2.gro -o no_h.gro -n no_h.ndx'
        p2 = subprocess.Popen(edit_cmd, shell=True, stdin=subprocess.PIPE,
                              stdout=output_stream, stderr=output_stream)
        p2.communicate('2\n')

    def _re_order(self, output_stream):
        # create a new index file
        lines = open('index.ndx').read().splitlines()
        header = lines[0]
        indices = []
        for line in lines[1:]:
            cols = line.split()
            for col in cols:
                indices.append( int(col) )
        resorted = [ indices.index(val)+1 for val in range( 1, max(indices)+1 ) ]
        with open('resort.ndx', 'w') as out:
            print >>out, header
            for val in resorted:
                print >>out, val

        # resort
        edit_cmd = join(GMX_PATH, 'editconf')
        edit_cmd += ' -f no_h.gro -o min.pdb -n resort.ndx'
        subprocess.check_call(edit_cmd, shell=True, stdout=output_stream,
                              stderr=output_stream)

    def run_minimization(self, posres_force_const=1000., output_stream=None):
        start = self._get_closest_frame()

        # create temp dir
        if os.path.isdir('Temp'):
            pass
        else:
            os.mkdir('Temp')
        os.chdir('Temp')

        # write the average file
        prody.writePDB('average.pdb', self.input_pdb)
        pdb_cmd = join(GMX_PATH, 'pdb2gmx')
        pdb_cmd += ' -f average.pdb -ff amber99sb-ildn -water none -n index.ndx -posrefc {} -o ref.gro -his'.format(
            posres_force_const)
        p = subprocess.Popen(pdb_cmd, shell=True, stdin=subprocess.PIPE,
                             stdout=output_stream, stderr=output_stream)
        p.communicate('0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n')
        # put it in a bigger box
        box_cmd = join(GMX_PATH, 'editconf')
        box_cmd += ' -f ref.gro -o ref_box.gro -c -box 999 999 999'
        subprocess.check_call(box_cmd, shell=True, stdout=output_stream,
                              stderr=output_stream)

        # write pdb file
        prody.writePDB('start.pdb', start)

        # pdb2gmx
        pdb_cmd = join(GMX_PATH, 'pdb2gmx')
        pdb_cmd += ' -f start.pdb -ff amber99sb-ildn -water none -n index.ndx -posrefc {} -his'.format(
            posres_force_const)
        p = subprocess.Popen(pdb_cmd, shell=True, stdin=subprocess.PIPE,
                             stdout=output_stream, stderr=output_stream)
        p.communicate('0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n0\n')

        # put it in a bigger box
        box_cmd = join(GMX_PATH, 'editconf')
        box_cmd += ' -f conf.gro -o box.gro -c -box 999 999 999'
        subprocess.check_call(box_cmd, shell=True, stdout=output_stream,
                              stderr=output_stream)

        #
        # Round 1
        #

        # write mdp file
        with open('min_round_1.mdp', 'w') as min_file:
            min_file.write( mdp_string.format(integrator='steep') )

        # run grompp
        grompp_cmd = join(GMX_PATH, 'grompp')
        grompp_cmd += ' -f min_round_1.mdp -c box.gro -p topol.top -o min_round_1 -r ref_box.gro'
        subprocess.check_call(grompp_cmd, shell=True, stdout=output_stream,
                              stderr=output_stream)

        # run mdrun
        md_cmd = join(GMX_PATH, 'mdrun')
        md_cmd += ' -deffnm min_round_1 -v -nt 1'
        subprocess.check_call(md_cmd, shell=True, stdout=output_stream,
                              stderr=output_stream)

        #
        # Round 2
        #

        # write mdp file
        with open('min_round_2.mdp', 'w') as min_file:
            min_file.write( mdp_string.format(integrator='l-bfgs') )

        # run grompp
        grompp_cmd = join(GMX_PATH, 'grompp')
        grompp_cmd += ' -f min_round_2.mdp -c min_round_1.gro -p topol.top -o min_round_2 -maxwarn 1 -r ref_box.gro'
        subprocess.check_call(grompp_cmd, shell=True, stdout=output_stream,
                              stderr=output_stream)

        # run mdrun
        md_cmd = join(GMX_PATH, 'mdrun')
        md_cmd += ' -deffnm min_round_2 -v -nt 1'
        subprocess.check_call(md_cmd, shell=True, stdout=output_stream,
                              stderr=output_stream)

        #
        # gather results
        #
        self._create_no_h_file(output_stream)
        self._re_order(output_stream)

        # load the pdb
        protein = prody.parsePDB('min.pdb').select('not hydrogen')

        # clean up
        os.chdir('..')
        shutil.rmtree('Temp')

        return protein


def main():
    r = parse_args()
    input_pdb_filename = r[0]
    trajectory_filename = r[1]
    output_pdb_filename = r[2]
    force_const = r[3]
    m = Minimizer(input_pdb_filename, trajectory_filename)
    minimized_protein = m.run_minimization(force_const)
    prody.writePDB(output_pdb_filename, minimized_protein)

if __name__ == '__main__':
    main()

