"""
Calculates salt bridges or electrostatic between two chains.

Arguments:
    saltBridges - use it to get stats for pair of 2 charges in distance of less than 4A.
    periphery - use it to look for periphery charges only

See also:
DRIEDING: A Generic Force field for molecular simulations
"""
from __future__ import print_function
import math

from . import DBConfig
from .kdtree import KDTree
from .ASA import ASA

VERBOSE = True

def getHbonds(pdb, pdbName):
    conn = DBConfig.get_connection()
    cursor = conn.cursor()
    cursor.execute("""select DonorChain,DonorResId,DonorSymbol,AccChain,AccResId,AccSymbol from
                        Ndrieding
                        inner join donors2
                        on DonorSymbol=donors2.Symbol
                        where PDB='%s'""" % pdbName)

    hbonds = set()
    for dChain, dResId, dSymbol, aChain, aResId, aSymbol in cursor.fetchall():
        donorAtom, accAtom = None, None
        for a in pdb.atoms:
            if a.chain == dChain and a.resId == dResId and a.symbol == dSymbol:
                donorRes = a.resId
            elif a.chain == aChain and a.resId == aResId and a.symbol == aSymbol:
                accRes = a.resId
        hbonds.add((donorRes, accRes))
    return hbonds

def fullHbonds(hbonds,charged_atom):
    charged_atomHbonds= [Hbond for Hbond in hbonds if charged_atom.resId in Hbond]
    if (charged_atom.residue in ['ASP', 'GLU'] or charged_atom.symbol == 'OXT') and len(charged_atomHbonds) >= 5:
        return fullHbonds
    if charged_atom.residue == 'LYS' and len(charged_atomHbonds) >= 3:
        return fullHbonds
    if charged_atom.residue == 'ARG' and len(charged_atomHbonds) >= 6:
        return fullHbonds


def eInteraction(Qi, Qj, R):
    kcal_mol_constant = 322.0637
    return kcal_mol_constant * Qi * Qj / (R ** 2)


def assignCharge(atom, pH=7):
    # how do we assign?

    if atom.residue in ['ASP', 'GLU'] and atom.atomType == 'O' and atom.symbol != 'O':
        return -0.5  # -1
    if atom.symbol == 'OXT':
        return -0.5  # -1

    # arg deloclalized
    # ARG - NH1 NH2 (not NE and N)
    # LYS NZ
    # HIS ND1 NE2
    if atom.residue == 'LYS' and atom.symbol == 'NZ':
        return 1.0
    posRes = ['ARG']
    if 0.1 < pH <= 6:
        posRes.append('HIS')
    if atom.residue in posRes and atom.atomType == 'N' and atom.symbol not in ['N', 'NE']:
        return 0.5  # 1

    return 0


def calcElectrostatic(pdb, interface, exclude_hbonds=False, count_cutoff_distance=5):
    """
    Calculated possible electro interactions, excluding 1-2 and 1-3 interactions
    (already included in angle and bond interactions
     """
    CUTOFF_DISTANCE = 7  # we could have 6?
    if exclude_hbonds:
        hHbonds = getHbonds(pdb, pdb.name)
    else:
        hHbonds = []

    components = []
    oxtAtoms = [(a.chain, a.resId) for a in pdb.atoms if a.symbol == 'OXT']  # C ter
    for a in pdb.atoms:
        if a.symbol == 'O' and (a.chain, a.resId) in oxtAtoms:
            a.symbol = 'OXT'

    for part in pdb.interfaceParts:
        components.append([a for a in interface if a.chain in part and assignCharge(a) != 0])

    # fill kd tree with charged atoms from second component
    comp2 = [atom for atom in pdb.atoms if atom.chain in pdb.interfaceParts[1] and assignCharge(atom) != 0]

    ktree = KDTree.construct_from_data(comp2)
    electroStat = 0.0
    pp, pm, mm = 0, 0, 0
    for atom in components[0]:
        # interactions are not calculated between atoms bonded to each other (1,2 and 1,3 [hbonds])
        Qi = assignCharge(atom)
        if Qi == 0:
            continue

        nearAtoms = list(ktree.findByDistance(query_point=atom.coord, distance=CUTOFF_DISTANCE ** 2))
        contact = [con for con in nearAtoms if not (((con, atom) in hHbonds) or ((atom, con) in hHbonds))]
        for con in contact:
            Qj = assignCharge(con)
            if Qj == 0:
                continue
            R = math.sqrt(atom.distance(con))
            electroStat += eInteraction(Qi, Qj, R)

            if R < count_cutoff_distance:
                if Qi > 0 and Qj > 0:
                    pp += 1
                elif Qi < 0 and Qj < 0:
                    mm += 1
                else:
                    pm += 1

    return electroStat, pp, mm, pm


def is_hydrophilic(atom):
    """
    Checks whether an atom belongs to hydrophilic residue
    :param atom: atom
    :return: True if the atom belongs to hydrophobic residue, otherwise false
    """

    hydrophilic_atoms = ['H', 'N', 'S', 'O']
    hydrophilic_residues = ['GLU', 'ASP', 'ASN', 'QLN', 'HIS', 'GLN', 'SER']
    
    if atom.symbol in ['CD', 'CZ'] and atom.residue == 'ARG':
        return is_hydrophilic
    if atom.symbol in ['CE2', 'CD1'] and atom.residue == 'TRP':
        return is_hydrophilic    
    if atom.symbol == 'CD' and atom.residue == 'PRO':
        return is_hydrophilic
    if atom.symbol == 'CB' and atom.residue == 'MET':
        return is_hydrophilic
    if atom.symbol in ['CD', 'CE'] and atom.residue == 'LYS':
        return is_hydrophilic
    if atom.symbol == 'CZ' and atom.residue == 'TYR':
        return is_hydrophilic
    if atom.symbol != 'CE' and atom.residue == 'MET':
        return is_hydrophilic
    if atom.symbol == 'C' or atom.symbol == 'CA':
        return is_hydrophilic
    if atom.residue in hydrophilic_residues:
        return is_hydrophilic
    if atom.atomType in hydrophilic_atoms:
        return is_hydrophilic

def charged_atom_areas(atom, pdb):
    
    if atom.symbol == 'OXT':
        charged = [a for a in pdb.atoms if a.symbol == 'OXT' and a.resId == atom.resId]
    if atom.residue == 'ASP':
        charged = [a for a in pdb.atoms if a.symbol in ['OD1', 'OD2'] and a.resId == atom.resId]
    if atom.residue == 'GLU':
        charged = [a for a in pdb.atoms if a.symbol in ['OE1', 'OE2'] and a.resId == atom.resId]
    if atom.residue == 'LYS':
        charged = [a for a in pdb.atoms if a.symbol in ['1HZ', '2HZ', '3HZ', 'NZ'] and a.resId == atom.resId]
    if atom.residue == 'ARG':
        charged = [a for a in pdb.atoms if a.symbol in ['NH1', 'NH2', 'NE'] and a.resId == atom.resId]
    return charged

def get_residue_depth(atom, depthDistances):
    for a, dist in depthDistances:
        if a.resId == atom.resId and a.symbol == atom.symbol:
            return dist

def calcInterElectroHydrophobic(pdb, interface, depthDistances):
    """
    Calculated possible electro interactions, excluding 1-2 and 1-3 interactions
    (already included in angle and bond interactions
     """
    HYDROPHOBIC_CHARGED_CUTOFF_DISTANCE = 6  # we could have 6?
    
    Hbonds = getHbonds(pdb, pdb.name)
    
    oxtAtoms = [(a.chain, a.resId) for a in pdb.atoms if a.symbol == 'OXT']  # C ter
    for a in pdb.atoms:
        if a.symbol == 'O' and (a.chain, a.resId) in oxtAtoms:
            a.symbol = 'OXT'

    hydrophobic_charged_interactions = []
    hydroElectroOutput = pdb.getFile('.hydrophobicElectro.txt') if VERBOSE else None

    for part in pdb.interfaceParts:
        charged_atoms = [a for a in interface if assignCharge(a) != 0 and a.chain in part and get_residue_depth(a, depthDistances) > 2.0]
        if not any(charged_atoms):
            continue
        electro_kdtree = KDTree.construct_from_data(charged_atoms)
        other_parts = ''.join([partb for partb in pdb.interfaceParts if partb != part])
        hydrophobic_partners = [a for a in interface if a.chain in other_parts and not is_hydrophilic(a)]
        for atom in hydrophobic_partners:
            nearAtoms = list(electro_kdtree.findByDistance(query_point=atom.coord, distance=HYDROPHOBIC_CHARGED_CUTOFF_DISTANCE ** 2))
            for con in nearAtoms:
                Qi = assignCharge(con)
                R = math.sqrt(atom.distance(con))
                depth = get_residue_depth (con, depthDistances)
                if R < HYDROPHOBIC_CHARGED_CUTOFF_DISTANCE and not (((con.resId, atom.resId) in Hbonds) or ((atom.resId, con.resId) in Hbonds)):
                    if Qi > 0:
                        hydrophobic_charged_interactions.append((con, depth, 'positive', atom))
                    elif Qi < 0:
                        hydrophobic_charged_interactions.append((con, depth, 'negative', atom))
                    if hydroElectroOutput:
                        print(','.join((con.chain, str(con.resId), con.residue, atom.chain, str(atom.resId),
                                       atom.residue, '%.3f' % R)), file=hydroElectroOutput)
    if hydroElectroOutput:
        hydroElectroOutput.close()
    return hydrophobic_charged_interactions    

def calcIntraElectroHydrophobic(pdb, interface, depthDistances):
    """
    Calculated possible electro interactions, excluding 1-2 and 1-3 interactions
    (already included in angle and bond interactions
     """
    HYDROPHOBIC_CHARGED_CUTOFF_DISTANCE = 6  # we could have 6?
    
    Hbonds = getHbonds(pdb, pdb.name)
    
    oxtAtoms = [(a.chain, a.resId) for a in pdb.atoms if a.symbol == 'OXT']  # C ter
    for a in pdb.atoms:
        if a.symbol == 'O' and (a.chain, a.resId) in oxtAtoms:
            a.symbol = 'OXT'

    hydrophobic_charged_interactions = []
    hydroElectroOutput = pdb.getFile('.hydrophobicElectro.txt') if VERBOSE else None

    for part in pdb.interfaceParts:
        charged_atoms = [a for a in interface if assignCharge(a) != 0 and a.chain in part and
                         get_residue_depth(a, depthDistances) > 2.0]
        if not any(charged_atoms):
            continue
        electro_kdtree = KDTree.construct_from_data(charged_atoms)
        hydrophobic_partners = [a for a in interface if a.chain in part and not is_hydrophilic(a)]
        for atom in hydrophobic_partners:
            nearAtoms = list(electro_kdtree.findByDistance(query_point=atom.coord, distance=HYDROPHOBIC_CHARGED_CUTOFF_DISTANCE ** 2))
            for con in nearAtoms:
                depth = get_residue_depth(con, depthDistances)
                Qi = assignCharge(con)
                R = math.sqrt(atom.distance(con))
                if R < HYDROPHOBIC_CHARGED_CUTOFF_DISTANCE and atom.resId != con.resId and not (((con.resId, atom.resId) in Hbonds) or ((atom.resId, con.resId) in Hbonds)):
                    if Qi > 0:
                        hydrophobic_charged_interactions.append((con, depth, 'positive', atom))
                    elif Qi < 0:
                        hydrophobic_charged_interactions.append((con, depth, 'positive', atom))
                    if hydroElectroOutput:
                        print(','.join((con.chain, str(con.resId), con.residue, atom.chain, str(atom.resId),
                                       atom.residue, '%.3f' % R)), file=hydroElectroOutput)
    if hydroElectroOutput:
        hydroElectroOutput.close()
    return hydrophobic_charged_interactions    
