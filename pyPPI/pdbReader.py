"""
Class for reading PDB files

References:
*http://www.wwpdb.org/documentation/format33/sect9.html#ATOM -pdb RFC
"""

import logging
import math
import os

from .atom import atom, water
from .kdtree import KDTree

NEIGHBOR_DISTANCE = 16
PDBS_DIR = None

def angle(DH, distance, HA, radians=True):
    ang = math.acos((DH ** 2 - distance ** 2 + HA ** 2) / (2 * DH * HA))
    if not radians:
        ang *= 180 / math.pi
    return ang


def radianToAngle(radian):
    """Get angle from radian
    :param radian: Radian
    :return: angle
    """
    return radian * (180 / math.pi)


class PDBReader(object):
    """ class that handles PDB files and utilities """

    @staticmethod
    def readFile(path, interface_parts=None):
        """Reads PDB file
        :param path:  path of the file
        :param interface_parts: relevant chains
        :return: a PDB object
        """
        name = os.path.basename(path)[0:4]
        pdbFile = open(path, 'r')
        atoms = []
        waters = []

        cmpnds = []
        cmpndChains = []
        hetAtms = []
        logging.info('reading pdb file (atoms and HETATM of HOH) %s', path)
        # read file for chain atoms and water
        for line in pdbFile.readlines():
            if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':
                atomNum, atomSymbol = line[6:11].strip(), line[12:16].strip()
                residue, chain, resId = line[17:20].strip(), line[21], line[22:26].strip(' '),
                x, y, z, atomType = line[30:38], line[38:46], line[46:54], line[77:78]
                tempFactor = float(line[60:66])

                # skip: when two rotamers appear, such as AVAL and BVAL we ignore the BVAL
                if line[16:17] != ' ' and line[16:17] != 'A':
                    continue
                # skip on multimodel. see 1EAW Res 60 chain A
                if line[26:27] != ' ':
                    continue
                coord = (float(x), float(y), float(z))
                if line[0:6] == 'ATOM  ':
                    atoms.append(atom(atomNum, atomSymbol, residue, chain, resId, atomType, coord, beta_factor=tempFactor))
                elif line[0:6] == 'HETATM':
                    if residue == 'HOH':
                        waters.append(water(atomNum, atomSymbol, residue, chain, resId, atomType, coord))
                    else:
                        hetAtms.append(atom(atomNum, atomSymbol, residue, chain, resId, atomType, coord, beta_factor=tempFactor))
            elif line[0:6] == "COMPND" and 'MOLECULE:' in line:
                cmpnds.append(line.split('MOLECULE:')[1].strip('; \n').replace(',', ' '))
            elif line[0:6] == "COMPND" and 'CHAIN:' in line:
                cmpndChains.append(line.split('CHAIN: ')[1].strip('; \n').replace(', ', ''))
            elif line[0:14] == "REMARK 200  PH":
                pH = line[45:50]
            elif line[0:3] == "END" and len(line) == 4:
                break
            elif line[0:14] == "MODEL        2":
                logging.info('using model 1 (more models are ignored)')
                break
        logging.info('finished reading file. waters:' + str(len(waters)) + ' atoms: ' + str(len(atoms)))
        if interface_parts is None:
            interface_parts = cmpndChains
        return PDBReader(name, atoms, waters, cmpnds, hetAtms, interface_parts=interface_parts)

    def __init__(self, name, atoms, waters, compunds, hetAtms, interface_parts=['A', 'B']):
        self.compunds = ' - '.join(compunds)
        self.interfaceParts = interface_parts
        self.name = name
        self.atoms = atoms
        self.waters = waters
        self.hetatms = hetAtms
        self.chains = set([a.chain for a in atoms])
        self.interfaceCache = None
        self.cacheDistance = 0
        self.__buildIndex()

    def __buildIndex(self):
        """ Init the internal indexs for the atoms and their pseudoChain, and k-tree """
        logging.debug('building k-tree')
        self.ktree = KDTree.construct_from_data(self.atoms[:])
        logging.debug('end building k-tree')
        logging.debug('building indexs')
        for i, a in enumerate(self.atoms):
            a.atomIndex = i
            # assign pseudo chain to residue which symbols the interface chains (example: A:B)
            for j, interPart in enumerate(self.interfaceParts):
                if a.chain in interPart:
                    a.pseudoChain = chr(ord('A') + j)

        logging.debug('end building index')

    def getFile(self, name):
        if PDBS_DIR is None:
            path = "./debug/"
        else:
            path = os.path.join(PDBS_DIR, 'debug')
        if not os.path.exists(path):
            os.makedirs(path)

        return open(os.path.join(path, self.name + name), 'w')

    def getPseudoChains(self):
        for j, interPart in enumerate(self.interfaceParts):
            yield chr(ord('A') + j)

    def getInterface(self, max_distance=NEIGHBOR_DISTANCE):
        """Get atoms not from same chain, having distance less than maxDistance ignores H
        """
        if self.interfaceCache is not None and self.cacheDistance == max_distance:
            return self.interfaceCache

        self.cacheDistance = max_distance
        self.interfaceCache = set()
        interfacesT = self.interfaceCache
        # we assume dimers

        for atom in (atom for atom in self.atoms if atom.pseudoChain == 'A'):
            for atom2 in self.ktree.findByDistance(query_point=atom.coord, distance=max_distance):
                if atom2.pseudoChain == atom.pseudoChain or atom2.pseudoChain is None:
                    continue
                interfacesT.add(atom)
                interfacesT.add(atom2)

        # extend with H
        print('interface atoms:', len(interfacesT))
        return self.interfaceCache

    def atoms(self):
        """Get atoms in the PDB"""
        return self.atoms

    def getNextAtoms(self, atom, i):
        """Get the next atom

        :param atom: atom to start from
        :param i: offset from the atom
        :return: the i'th atom from atom
        """
        if atom.atomIndex + i >= len(self.atoms):
            logging.debug('error getting atom %s', atom)
            return atom
        return self.atoms[atom.atomIndex + i]
