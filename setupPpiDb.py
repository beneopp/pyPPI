#!/usr/bin/env python
"""
Setups a protein database in MySQL: a database of interesting properties of the proteins based on scripts of this library.

This should be easy to use script for invoking the most important scripts of the library and store them in DB
for easy retrieve.

How to use:
Create a folder and place there some file with list of PDBs to analyze.
The program will create the following directory structure in the same directory:
    ./pdbs/ - list of pdbs downloaded
    ./results/ - results of the analysis scripts
"""
from __future__ import print_function
import argparse
import os
import subprocess
import sys
import pkg_resources
import requests
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.ResidueDepth import get_surface, min_dist

from pyPPI import DBConfig
import pyPPI.surfaceComplementarity.VDW as VDW
import pyPPI.surfaceComplementarity.interfaceDepth as Periphery
from pyPPI.ASA import ASA
from pyPPI.hbonds import hbonds
from pyPPI.kdtree import KDTree
import pyPPI.pdbReader as pdbReader
from pyPPI.pdbReader import PDBReader
import pyPPI.electrostat as electrostat
from pyPPI.cavities import calculateVolume


"""
Distance in angtroms between the chains that is relevant for defining the interface
"""
INTERFACE_DISTANCE = 4
WORKING_DIRECTORY = './'
PDBS_DIR = "./pdbs/"
RESULTS_DIR = "./results/"

_remediator = pkg_resources.resource_filename('pyPPI', '/'.join(['molprobity', 'remediator.pl']))
_reduce_path = pkg_resources.resource_filename('pyPPI', '/'.join(['molprobity', 'reduce']))

def download_PDB(pdb):
    """
    Downloads a PDB from protein data base
    :param pdb: pdb identifier
    """
    url = 'http://www.rcsb.org/pdb/files/{0}.pdb'.format(pdb)
    print('downloading %s (%s)' % (pdb, url))

    req = requests.get(url)
    with get_file(pdb) as newPDB:
        print(req.text, file=newPDB)


def get_file(name):
    """
    Get file for write in the PDBS_DIR
    :param name:
    :return:
    """
    global PDBS_DIR
    return open(os.path.join(PDBS_DIR, name + ".pdb"), "w")


def download_DB(pdbList):
    """
    Downloads PDB and add hydrogens using molprobity
    :param pdbList: list of pdbs to download
    """
    print("Downloading pdbs according to list")
    for pdb in pdbList:
        # don't download twice the same PDB
        if os.path.exists(os.path.join(PDBS_DIR, pdb + "_FH.pdb")): continue

        # in case the PDB is already in the directory
        if not os.path.exists(os.path.join(PDBS_DIR, pdb + ".pdb")):
            download_PDB(pdb)

        molprobity(pdb)
    print("Finished downloading pdbs")


def molprobity(pdb_name):
    """
    runs molprobility on a input protein
    :param pdb_name: name of the PDB file
    :return:
    """
    global MOLPROBITY_DIR, PDBS_DIR
    if os.path.exists(os.path.join(PDBS_DIR, pdb_name + "_FH.pdb")):
        return True  # already exist
    print('Starting molprobity %s' % pdb_name)
    subprocess.check_output('perl ' + _remediator + ' ' + os.path.join(PDBS_DIR,
                                                                                                         pdb_name + ".pdb") + ' > a',
                            shell=True)
    try:
        subprocess.check_output(_reduce_path + ' a > b', shell=True)
    except:
        print('error prasing PDB %s' % pdb_name)
        pass  # yakky kaky, but reduce returns 1 exit
    subprocess.check_output(
        'perl ' + _remediator +' b -oldout> ' + os.path.join(PDBS_DIR, pdb_name + "_FH.pdb"),
        shell=True)
    # delete the PDB file - we will work with a file with hydrogens added (_FH create above)
    os.remove(os.path.join(PDBS_DIR, pdb_name + ".pdb"))


def buildASAperAtomForComplex(pdb, result):
    asaCalc = ASA(pdb)
    asaCalc.execute()
    for atom, asa in asaCalc.interPerAtom.items():
        # complex inter
        res = [pdb.name, atom.chain, atom.residue, atom.resId, atom.symbol, atom.atomType, asa, atom.tempFactor, 0]
        print(','.join([str(a) for a in res]), file=result)
        # complex intra (separated)
        asa = asaCalc.diffASAperAtom[atom] + asa
        res = [pdb.name, atom.chain, atom.residue, atom.resId, atom.symbol, atom.atomType, asa, atom.tempFactor, 1]
        print(','.join([str(a) for a in res]), file=result)

def calcInterfaceDist(pdb, result):
    """
        Defines interface by distance
        """
    global INTERFACE_DISTANCE
    partA = [a for a in pdb.atoms if a.chain in pdb.interfaceParts[0]]
    partB = [a for a in pdb.atoms if a.chain in pdb.interfaceParts[1]]
    if len(partA) == 0 or len(partB) == 0:
        print('WARNING: %s doesnt have atoms in one its chains' % pdb.name)
        return
    aTree = KDTree.construct_from_data(partA[:])
    bTree = KDTree.construct_from_data(partB[:])
    complexChains = ':'.join(pdb.interfaceParts)
    for part, tree in [(partA, bTree), (partB, aTree)]:
        for atom in part:
            near, dist = tree.findNearest(query_point=atom.coord, num=1)
            if dist < INTERFACE_DISTANCE:
                print(','.join([pdb.name, complexChains, atom.chain, str(atom.resId), atom.symbol, atom.atomType, str(dist)]), file=result)

def createInterfaceCSV(pdbsToAnalyze):
    """
        interface can be defined by either ASA or distance
        we use both of them
        """
    global PDBS_DIR, RESULTS_DIR
    if all(os.path.exists(os.path.join(RESULTS_DIR, resFile)) for resFile in ['PerAtomASA.csv', 'PerAtomASA.csv']):
        print('Data already exist in result directory.')
        return
    
    with open(os.path.join(RESULTS_DIR, 'PerAtomASA.csv'), 'w') as asaPerAtom:
        with open(os.path.join(RESULTS_DIR, 'PerAtomDistance.csv'), 'w') as distancePerAtom:
            pdbs = os.listdir(PDBS_DIR)
            print('PDB,Chains,Chain,ResId,Symbol,Atom,MinDistance', file=distancePerAtom)
            print('PDB,Chain,Residue,ResId,Symbol,AtomType,ASA,tempFactor,Seperated', file=asaPerAtom)
            failedPDBs = []
            pdbsNamesToChains = dict((p[0], p[1].split(':') if len(p) > 1 else None) for p in pdbsToAnalyze)
            for pdbName in pdbs:
                if pdbName[0:4] not in pdbsNamesToChains: continue
                pdb = PDBReader.readFile(os.path.join(PDBS_DIR, pdbName), pdbsNamesToChains[pdbName[0:4]])
                try:
                    print('Writing ASA for %s' % pdb.name)
                    buildASAperAtomForComplex(pdb, asaPerAtom)
                    print('Writing distance for %s' % pdb.name)
                    calcInterfaceDist(pdb, distancePerAtom)
                except IndexError:
                    failedPDBs.append(pdb.name)

    print('Finished')
    if len(failedPDBs) > 0:
        print('Failed to process:', ','.join(failedPDBs))


def createDataBase(pdbsToAnalyzeWithChains):
    """Loads teh computations to a new database
    :param pdbsToAnalyzeWithChains:
    """
    print('Creating DB: %s' % DBConfig.DB_NAME)

    installDB = pkg_resources.resource_filename('pyPPI', '/'.join(['sqls', 'createDB.sql']))
    metadataDB = pkg_resources.resource_filename('pyPPI', '/'.join(['sqls', 'donors2.sql']))
    createInterfaceSql = pkg_resources.resource_filename('pyPPI', '/'.join(['sqls', 'createInterface.sql']))

    subprocess.call(
        "mysql -u%s -p%s -e 'create database if not exists %s'" % (DBConfig.USER, DBConfig.PASSWD, DBConfig.DB_NAME),
        shell=True)
    # create schema
    subprocess.call('mysql %s -u%s -p%s < %s ' % (DBConfig.DB_NAME, DBConfig.USER, DBConfig.PASSWD, installDB),
                    shell=True)
    # insert metadata
    subprocess.call('mysql %s -u%s -p%s < %s ' % (DBConfig.DB_NAME, DBConfig.USER, DBConfig.PASSWD, metadataDB),
                    shell=True)
    conn = DBConfig.get_connection()
    cursor = conn.cursor()
    cursor.execute('''
    load data local infile '%s' into table interfaceDist fields terminated by ',' optionally enclosed by '"'  lines terminated by '\n' ignore 1 lines (PDB,Chains,Chain,ResId,Symbol,Atom,MinDist);
    ''' % (os.path.join(RESULTS_DIR, 'PerAtomDistance.csv')))
    cursor.execute('''
    load data local infile '%s' into table perAtomASA fields terminated by ',' optionally enclosed by '"'  lines terminated by '\n' ignore 1 lines (PDB,Chain,Residue,ResId,Symbol,Atom,ASA,Bfactor,Seperated);
    ''' % (os.path.join(RESULTS_DIR, 'PerAtomASA.csv')))
    conn.commit()

    # create interface table
    subprocess.call('mysql %s -u%s -p%s < %s ' % (DBConfig.DB_NAME, DBConfig.USER, DBConfig.PASSWD, createInterfaceSql),
                    shell=True)

    # add metadata table with complexs in the database
    pdbsNamesToChains = dict((p[0], p[1].split(':') if len(p) > 1 else None) for p in pdbsToAnalyzeWithChains)
    dataToInsert = []
    for pdbName, chains in pdbsNamesToChains.items():
        pdb = PDBReader.readFile(os.path.join(PDBS_DIR, '%s_FH.pdb' % pdbName), pdbsNamesToChains[pdbName[0:4]])
        if chains is None:
            compunds = pdb.compunds.split(' - ')
            dataToInsert.append((pdbName, pdb.interfaceParts[0], compunds[0] if len(compunds) > 1 else compunds,
                                 pdb.interfaceParts[1], compunds[1] if len(compunds) > 1 else ''))
        else:
            dataToInsert.append((pdbName, pdb.interfaceParts[0], '', pdb.interfaceParts[1], ''))

    cursor = conn.cursor()
    cursor.executemany('''
    INSERT INTO proteinComplex (PDB,UnboundChainA,NameA,UnboundChainB,NameB)
    values (%s,%s,%s,%s,%s)
    ''', dataToInsert)
    conn.commit()
    conn.close()
    print('database created!')


def getInterfaceAtoms(cur, pdb):
    """
    Gets interface atoms from database
    :param cur: cursor to database
    :param pdb: pdb object to get atoms from
    :return: list of interface atoms
    """
    cur.execute('''
    select Chain,ResId,Symbol from NinterfaceAtoms
    where PDB='%s'
    ''' % pdb.name)
    interfaceAtoms = []
    for chain, resid, symbol in cur.fetchall():
        interfaceAtoms.append(
            next(a for a in pdb.atoms if a.chain == chain and a.resId == resid and a.symbol == symbol))
    return interfaceAtoms


def fillInterfacePeriphrial(pdbsToAnalyze):
    global PDBS_DIR, RESULTS_DIR

    if os.path.exists(os.path.join(RESULTS_DIR, 'interfacePeriphrial.csv')):
        print('Data already exist in result directory for interface periphery.')
        return

    pdbsNamesToChains = dict((p[0], p[1].split(':') if len(p) > 1 else None) for p in pdbsToAnalyze)

    with open(os.path.join(RESULTS_DIR, 'interfacePeriphrial.csv'), 'w') as interfacePeriphrial:
        print('PDB,Chain,ResId,Symbol,Peripherial,PropPeri', file=interfacePeriphrial)
        for pdbName, chains in pdbsNamesToChains.items():
            print('Calculating peripheral table for %s ' % pdbName)
            pdb_path = os.path.join(PDBS_DIR, '%s_FH.pdb' % pdbName)
            depthL, peripherialL = Periphery.calc_peripheral_PDB(pdb_path, chains)
            for atom, peri, propPeri in peripherialL:
                print(','.join([pdbName, atom.chain, str(atom.resId), atom.symbol, str(peri), str(propPeri)]),
                      file=interfacePeriphrial)

    conn = DBConfig.get_connection()
    cursor = conn.cursor()
    cursor.execute('''
    load data local infile '%s' into table interfacePeriphrial
    fields terminated by ',' optionally enclosed by '"'  lines terminated by '\n'
    ignore 1 lines (PDB,Chain,ResId,Symbol,Peri,PropPeri);
    ''' % (os.path.join(RESULTS_DIR, 'interfacePeriphrial.csv')))
    conn.commit()
    conn.close()

    conn = DBConfig.get_connection()
    cursor = conn.cursor()
    cursor.execute('''
    load data local infile '%s' into table interfacePeriphrial
    fields terminated by ',' optionally enclosed by '"'  lines terminated by '\n'
    ignore 1 lines (PDB,Chain,ResId,Symbol,Peri,PropPeri);
    ''' % (os.path.join(RESULTS_DIR, 'interfacePeriphrial.csv')))
    conn.commit()
    conn.close()

def residue_depth(pdbName, ReaderAtomsInput, filename, UseInterfaceAtoms=False):
    
    parser = PDBParser(PERMISSIVE=1)
    structure = parser.get_structure(pdbName, filename)
    model = structure[0]

    BioAtoms = []
    for chain in model:
        for residue in chain:
            for atom in residue:
                BioAtoms.append(atom)
    if UseInterfaceAtoms:
        BioAtoms = pdbReader_to_BioPyth(ReaderAtomsInput, BioAtoms)

    surface = get_surface(model)
    BioDepthDistances = []
    for atom in BioAtoms:
        dist = min_dist(atom.get_coord(), surface)
        BioDepthDistances.append([atom, dist])
    pdbReaderDistances = BioPyth_to_pdbReader(BioDepthDistances, ReaderAtomsInput)
    return pdbReaderDistances

def pdbReader_to_BioPyth(ReaderAtomsInput, BioAtomsInput):
    BioAtomsOutput = []
    
    for a in ReaderAtomsInput:
        for atom in BioAtomsInput:
            residue = atom.get_parent()
            if residue.get_id() == (' ', a.resId, ' ') and atom.get_name() == a.symbol:
                BioAtomsOutput.append(atom)
                break
    return BioAtomsOutput

def BioPyth_to_pdbReader(BioAtomsInput, ReaderAtomsInput):
    ReaderAtomsOutput = []
    for atom, dist in BioAtomsInput:
        residue = atom.get_parent()
        for a in ReaderAtomsInput:
            if residue.get_id() == (' ', a.resId, ' ') and atom.get_name() == a.symbol:
                ReaderAtomsOutput.append([a, dist])
                break
    return ReaderAtomsOutput

def calcEnergyTerms(pdbsToAnalyze):
    """
    Finds hydrogen bonds near interface atoms and calculates their energy,
    and calculates VDW and electrostatic energy for PDB
    """
    global PDBS_DIR, RESULTS_DIR
    
    output_file_list = ['Ndrieding.csv', 'interfaceVDW.csv', 'electrostatic.csv', 'electrostatic-hydrophobic.csv', 'cavity_vol.csv', 'cavity_res.csv', 'residue_depth.csv']
    if all(os.path.exists(os.path.join(RESULTS_DIR, resFile)) for resFile in output_file_list):
        print('Data already exists in result directory for energy terms.')
        return

    conn = DBConfig.get_connection()
    cursor = conn.cursor()
    pdbsNamesToChains = dict((p[0], p[1].split(':') if len(p) > 1 else None) for p in pdbsToAnalyze)

    if os.path.exists(os.path.join(RESULTS_DIR, 'Ndrieding.csv')) and os.path.exists(os.path.join(RESULTS_DIR, 'electrostatic-hydrophobic.csv')):
        print('Skipping hydrophobic calculations since file already exists.')
    if not os.path.exists(os.path.join(RESULTS_DIR, 'Ndrieding.csv')) or not os.path.exists(os.path.join(RESULTS_DIR, 'electrostatic-hydrophobic.csv')):
    #the Ndrieding.csv is needed if electrostatic-hydrophobic.csv does not exist because it loads the Hbonds into the database
        with open(os.path.join(RESULTS_DIR, 'Ndrieding.csv'), 'w') as driedingResult:
            print('PDB,DonorChain,DonorResId,DonorSymbol,AccChain,AccResId,AccSymbol,Energy', file=driedingResult)
            pdbs = os.listdir(PDBS_DIR)
            for pdbName in pdbs:
                if pdbName[0:4] not in pdbsNamesToChains: continue
                pdb = PDBReader.readFile(os.path.join(PDBS_DIR, pdbName), pdbsNamesToChains[pdbName[0:4]])
                interfaceAtoms = getInterfaceAtoms(cursor, pdb)
                bonds = hbonds(pdb)
                bonds.HDPlusDefinition = False
                cBondList = bonds.hbonds(interfaceAtoms)
                print('Calcing Hbonds for %s' % pdb.name)
                for donor, acceptor, eng in cBondList:
                    toPrint = [pdb.name, donor.chain, donor.resId, donor.symbol, acceptor.chain, acceptor.resId,
                               acceptor.symbol, eng]
                    print(','.join([str(a) for a in toPrint]), file=driedingResult)
        cursor.execute('''
        load data local infile '%s' into table Ndrieding
        fields terminated by ',' optionally enclosed by '"'  lines terminated by '\n'
        ignore 1 lines (PDB,DonorChain,DonorResId,DonorSymbol,AccChain,AccResId,AccSymbol,Energy);
        ''' % (os.path.join(RESULTS_DIR, 'Ndrieding.csv')))
        conn.commit()
                         
    if os.path.exists(os.path.join(RESULTS_DIR, 'interfaceVDW.csv')):
        print('Skipping VDW energy calculations since file already exists.')
    if not os.path.exists(os.path.join(RESULTS_DIR, 'interfaceVDW.csv')):
        print('Calculating VDW energy between interfaces')
        with open(os.path.join(RESULTS_DIR, 'interfaceVDW.csv'), 'w') as vdw_result:
            print('PDB,VDV,VDVx,clashV,clashS', file=vdw_result)
            for pdb, chains in pdbsNamesToChains.items():
                print('Calcing VDW for %s' % pdb)
                pdb_path = os.path.join(PDBS_DIR, '%s_FH.pdb' % pdb)
                sumVDW, sumVDWx, clashV, clashS = VDW.calcCompl(pdb_path, chains)
                print(','.join([pdb, str(sumVDW), str(sumVDWx), str(clashV), str(clashS)]), file=vdw_result)
        cursor.execute('''
        load data local infile '%s' into table interfaceVDW
        fields terminated by ',' optionally enclosed by '"'  lines terminated by '\n'
        ignore 1 lines (PDB,VDV,VDVx6,ClashV,ClashS);
        ''' % (os.path.join(RESULTS_DIR, 'interfaceVDW.csv')))
        conn.commit()
    
    if os.path.exists(os.path.join(RESULTS_DIR, 'electrostatic.csv')):
        print('Skipping electrostatic charge calculations since file already exists.')
    if not os.path.exists(os.path.join(RESULTS_DIR, 'electrostatic.csv')):
        print('Calculating electrostatic charges (Coulomb of paired charges except hydrogen bonds)')
        with open(os.path.join(RESULTS_DIR, 'electrostatic.csv'), 'w') as electro_res:
            print('PDB,eCoulomb,pp,mm,pm', file=electro_res)
            for pdb, chains in pdbsNamesToChains.items():
                pdb_path = os.path.join(PDBS_DIR, '%s_FH.pdb' % pdb)
                pdb = PDBReader.readFile(pdb_path, chains)
                interfaceAtoms = getInterfaceAtoms(cursor, pdb)

                e, pp, mm, pm = electrostat.calcElectrostatic(pdb, interfaceAtoms)
                print('%s,%f,%i,%i,%i' % (pdb.name, e, pp, mm, pm), file=electro_res)
        cursor.execute('''
        load data local infile '%s' into table electrostat
        fields terminated by ',' optionally enclosed by '"'  lines terminated by '\n'
        ignore 1 lines (PDB,electro,pp,mm,pm);
        ''' % (os.path.join(RESULTS_DIR, 'electrostatic.csv')))
        conn.commit()

    if os.path.exists(os.path.join(RESULTS_DIR, 'electrostatic-hydrophobic.csv')):
        print('Skipping electrostatic-hydrophobic interaction calculations since file already exists.')
    if not os.path.exists(os.path.join(RESULTS_DIR, 'electrostatic-hydrophobic.csv')):
        print('Calculating intermolecular electrostatic charges contacts with hydrophobic residues')
        with open(os.path.join(RESULTS_DIR, 'electrostatic-hydrophobic.csv'), 'w') as electro_hydro_res:
            print('inter-chain interactions', file=electro_hydro_res)
            print('\nPDB Name,charged atom chain,charged atom resID,charged atom residue,charged atom symbol,charged atom depth,charge,hydrophobic atom chain,hydrophobic atom resID,hydrophobic atom residue,hydrophobic atom symbol', file=electro_hydro_res)
            for pdbName, chains in pdbsNamesToChains.items():
                if pdbName[0:4] not in pdbsNamesToChains: continue
                pdb_path = os.path.join(PDBS_DIR, '%s_FH.pdb' % pdbName)
                pdb = PDBReader.readFile(pdb_path, chains)
                interfaceAtoms = getInterfaceAtoms(cursor, pdb)
                depthDistances = residue_depth(pdbName[:4], interfaceAtoms, pdb_path, UseInterfaceAtoms=True)
                inter_hydrophobic_charged_interactions = electrostat.calcInterElectroHydrophobic(pdb, interfaceAtoms, depthDistances)
                for charged_atom, depth, charge, hydrophobic_atom in inter_hydrophobic_charged_interactions:
                    print(','.join([pdbName, charged_atom.chain, str(charged_atom.resId), charged_atom.residue, charged_atom.symbol, str(depth), charge, hydrophobic_atom.chain, str(hydrophobic_atom.resId), hydrophobic_atom.residue, hydrophobic_atom.symbol]), file=electro_hydro_res)
    
            print('Calculating intramolecular electrostatic charges contacts with hydrophobic residues')
            print('\nintra-chain interactions', file=electro_hydro_res)
            print('\nPDB Name,charged atom chain,charged atom resID,charged atom residue,charged atom symbol,charged atom depth,charge,hydrophobic atom chain,hydrophobic atom resID,hydrophobic atom residue,hydrophobic atom symbol', file=electro_hydro_res)
            for pdbName, chains in pdbsNamesToChains.items():
                if pdbName[0:4] not in pdbsNamesToChains: continue
                pdb_path = os.path.join(PDBS_DIR, '%s_FH.pdb' % pdbName)
                pdb = PDBReader.readFile(pdb_path, chains)
                interfaceAtoms = getInterfaceAtoms(cursor, pdb)
                depthDistances = residue_depth(pdbName[:4], interfaceAtoms, pdb_path, UseInterfaceAtoms=True)
                intra_hydrophobic_charged_interactions = electrostat.calcIntraElectroHydrophobic(pdb, interfaceAtoms, depthDistances)
                for charged_atom, depth, charge, hydrophobic_atom in intra_hydrophobic_charged_interactions:
                    print(','.join([pdbName, charged_atom.chain, str(charged_atom.resId), charged_atom.residue, charged_atom.symbol, str(depth), charge, hydrophobic_atom.chain, str(hydrophobic_atom.resId), hydrophobic_atom.residue, hydrophobic_atom.symbol]), file=electro_hydro_res)
    
    if os.path.exists(os.path.join(RESULTS_DIR, 'cavity_vol.csv')):
        print('Skipping cavity calculations since file already exists.')
    if not os.path.exists(os.path.join(RESULTS_DIR, 'cavity_vol.csv')):
        print('Approximating cavities/gaps volume by monte carlo')
        with open(os.path.join(RESULTS_DIR, 'cavity_vol.csv'), 'w') as cavity_res:
            print('PDB,cavity_vol', file=cavity_res)
            pdbs = os.listdir(PDBS_DIR)
            for pdbName in pdbs:
                if pdbName[0:4] not in pdbsNamesToChains: continue
                pdb = PDBReader.readFile(os.path.join(PDBS_DIR, pdbName), pdbsNamesToChains[pdbName[0:4]])
                interfaceAtoms = getInterfaceAtoms(cursor, pdb)
                cavities_vol_approx = calculateVolume(pdb, interfaceAtoms)
                print('%s,%f' % (pdb.name, cavities_vol_approx), file=cavity_res)

    if os.path.exists(os.path.join(RESULTS_DIR, 'residue_depth.csv')):
        print('Skipping residue depth calculations since file already exists.')
    if not os.path.exists(os.path.join(RESULTS_DIR, 'residue_depth.csv')):
        print('Calculating Residue Depth')
        pdbs = os.listdir(PDBS_DIR)
        with open(os.path.join(RESULTS_DIR, 'residue_depth.csv'), 'w') as res_depth:
            print('PDB,chain,resId,residue,atom symbol,atom depth', file=res_depth)
            for pdbName in pdbs:
                if pdbName[0:4] not in pdbsNamesToChains: continue
                pdb_path = os.path.join(PDBS_DIR, '%s_FH.pdb' % pdbName[:4])
                pdb = PDBReader.readFile(pdb_path, pdbName[:4])
                depthDistances = residue_depth(pdbName[:4], pdb.atoms, pdb_path)
                for atom, dist in depthDistances:
                    print(','.join([pdbName[:4], atom.chain, str(atom.resId), atom.residue, atom.symbol, str(dist)]), file=res_depth)

    cursor.close()
    conn.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Setup/download protein database based on PDB")
    parser.add_argument("pdbList", help="A file with a list of PDB to download")
    parser.add_argument("--folder", help="Name of the folder to contain downloaded files")
    parser.add_argument("--dbName", help="Name of the database to create.")
    args = parser.parse_args()
    if args.pdbList is None:
        sys.exit("Please provide a file with list of PDBs to anaylze")

    WORKING_DIRECTORY = args.folder if args.folder is not None else os.path.dirname(os.path.abspath(args.pdbList))
    print('WORKING DIR: %s' % WORKING_DIRECTORY)

    PDBS_DIR = os.path.join(WORKING_DIRECTORY, 'pdbs')
    pdbReader.PDBS_DIR = PDBS_DIR
    RESULTS_DIR = os.path.join(WORKING_DIRECTORY, 'results')
    for dir in [PDBS_DIR, RESULTS_DIR]:
        if not os.path.exists(dir):
            os.mkdir(dir)

    pdbsToAnalyzeWithChains = [pdb.strip().upper().split("_") for pdb in open(args.pdbList, 'r') if
                           pdb[0:1] != '#']  # todo: add treatment for chains specificatin instad of [0:4]
    pdbsToAnalyze = [pdb[0] for pdb in pdbsToAnalyzeWithChains]
    download_DB(pdbsToAnalyze)  # download from PDB bank and add hydrogens
    createInterfaceCSV(pdbsToAnalyzeWithChains)  # define interface by distance and by asa
    print('''The script will now create DB. DB is required for extra calculations
        including VDW and hydrogen bonds
        ''')
    try:
        if args.dbName:
            DBConfig.DB_NAME = args.dbName
        DBConfig.init_connection()
        createDataBase(pdbsToAnalyzeWithChains)
        
        # post database creation scripts
        fillInterfacePeriphrial(pdbsToAnalyzeWithChains)
        calcEnergyTerms(pdbsToAnalyzeWithChains)

    except KeyboardInterrupt:
        print('DB will not be created. Use ./results table to see the results')
