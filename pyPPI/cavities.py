"""
Monte carlo estimation of cavity volume.
"""

import numpy as np
from scipy.linalg import norm
from scipy.spatial import cKDTree

from .ASA import radio_atom, KNOWN_RADIUS, R_WATER
from .DBConfig import get_connection


def getPeriperial(pdb, atoms):
    conn = get_connection()
    cursor = conn.cursor()

    # create interface table
    cursor.execute("""select Chain,ResId,Symbol,PropPeri from
                        interfacePeriphrial
                        where PDB='%s'""" % pdb.name)

    periDict = dict()
    for chain, resId, symbol, peri in cursor.fetchall():
        atom = [a for a in atoms if a.symbol == symbol and a.chain == chain and a.resId == resId][0]
        periDict[atom] = peri
    return periDict


def calculateVolume(pdb, interface):
    ACCURACY = 1.0 ** 3
    ACCURACY_FACTOR = 10.0
    components = [[a for a in interface if a.chain in part] for part in pdb.interfaceParts]

    partResIds = [set((a.chain, a.resId) for a in part) for part in components]
    interfaceRes = [[a for a in pdb.atoms if (a.chain, a.resId) in partResIds] for part, partResIds
                    in zip(pdb.interfaceParts, partResIds)]

    # find bounding box
    boundingBox = []
    for i in range(0, 3):
        axisCoords = [a.coord[i] for a in interface]
        boundingBox.append((min(axisCoords), max(axisCoords)))

    distances = [maxAxis - minAxis for minAxis, maxAxis in boundingBox]
    allVolume = np.prod(distances)
    pointsToCheck = ACCURACY_FACTOR * (allVolume / ACCURACY)

    aTree = cKDTree(np.array([a.coord for a in interfaceRes[0]]))
    bTree = cKDTree(np.array([a.coord for a in interfaceRes[1]]))

    MAX_WDW_RAD = max(KNOWN_RADIUS.values())

    rand_points = np.random.rand(2 * int(pointsToCheck), 3)
    for i, min_max in enumerate(boundingBox):
        rand_points[:, i] *= min_max[1] - min_max[0]
        rand_points[:, i] += min_max[0]

    aDistances, neighbors_a = aTree.query(rand_points)
    bDistances, neighbors_b = bTree.query(rand_points)

    CAVITIY_DISTANCE_CUTOFF = 6

    nearA = np.array([a for a in interfaceRes[0]])[neighbors_a]
    nearB = np.array([a for a in interfaceRes[1]])[neighbors_b]
    ab_Distance = norm(np.array([a.coord for a in nearA]) - np.array([a.coord for a in nearB]), axis=1)
    # only points within the interface
    selector = (aDistances < ab_Distance) & (bDistances < ab_Distance)
    rand_points = rand_points[selector]
    aDistances, neighbors_a = aDistances[selector], neighbors_a[selector]
    bDistances, neighbors_b = bDistances[selector], neighbors_b[selector]

    radiusA = np.array([radio_atom(a.atomType) for a in nearA[selector]])
    radiusB = np.array([radio_atom(a.atomType) for a in nearB[selector]])
    interface_points = ((aDistances < (radiusA + 2 * R_WATER)) & (bDistances < (radiusB + 2 * R_WATER))).sum()

    selector = np.maximum(aDistances, bDistances) < CAVITIY_DISTANCE_CUTOFF
    rand_points = rand_points[selector]

    # sufficient volume to accommodate a water molecule
    interfaceKdTree = cKDTree(np.array([a.coord for a in interfaceRes[0] + interfaceRes[1]]))
    distance_rand_points, _ = interfaceKdTree.query(rand_points)
    cavities_candidates = (distance_rand_points >= (MAX_WDW_RAD + R_WATER)).sum()

    cavitiesVolume = (cavities_candidates / pointsToCheck) * allVolume
    return cavitiesVolume
