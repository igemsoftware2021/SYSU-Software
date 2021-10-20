#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from Bio import Align, SeqIO
from Bio.PDB.PDBParser import PDBParser
import numpy as np
from scipy.optimize import leastsq
import json

def getAligned(srcPDB, tarPDB):
    """
    Align the enzyme with the entire protein chain.

    Parameters
    ----------
    srcPDB : string
        The filename of the PDB file recording the enzyme.
    tarPDB : string
        The filename of the PDB file recording the entire chain (E1-linker-CRY2-CRY2-linker-E2).

    Returns
    -------
    aligned: tuple
        alignment.aligned obtained from PairwiseAligner, used to generate information about position correspondence.

    """
    srcRecord = SeqIO.read(srcPDB, "pdb-atom")
    tarRecord = next(SeqIO.parse(tarPDB, 'pdb-atom'))

    aligner = Align.PairwiseAligner()
    alignments = aligner.align(srcRecord.seq, tarRecord.seq)
    alignment = alignments[0]
    aligned = alignment.aligned
    return aligned


def getCoordinates(pdb, positions):
    """
    Obtain the coordinate of each atom and the residue on which it is located from the information about position correspondence.

    Parameters
    ----------
    pdb : string
        The filename of the PDB file.
    positions : list
        The information about position correspondence.

    Returns
    -------
    coordinates : list
        The coordinate of each atom.
    matchedResidues : list
        For each atom, the residue where it is located.    

    """
    parser = PDBParser(PERMISSIVE = 1)
    structure = parser.get_structure('CRY2', pdb)
    
    position_index = 0
    residue_index = 0
    coordinates = []
    matchedResidues = []
    
    for model in structure.get_models():
        for chain in model.get_chains():
            for residue in chain.get_residues():
                if residue_index == positions[position_index]:
                    for atom in residue.get_atoms():
                        coordinates.append(atom.get_coord())
                        matchedResidues.append(residue)
                position_index += 1
                if position_index == len(positions):
                    return coordinates, matchedResidues
                residue_index += 1
            
    print('Error')
    return None


def alignedToPositions(aligned):
    """
    Convert alignment.aligned obtained from PairwiseAligner a list recording the index of the atoms which are on the chain and come from the enzyme.
    """
    positions = []
    for subseq in aligned[1]:
        positions.extend([x for x in range(subseq[0], subseq[1])])
    return positions


def getAfterCoordinates(srcPDB, tarPDB):
    """
    Find on the chain the atoms that come from the enzyme.

    Parameters
    ----------
    srcPDB : string
        The filename of the PDB file recording the enzyme.
    tarPDB : string
        The filename of the PDB file recording the entire chain (E1-linker-CRY2-CRY2-linker-E2).

    Returns
    -------
    TYPE
        DESCRIPTION.
    np.array(coordinates) : numpy Array of float64
        The coordinate of each atom from the enzyme.
    positions : list
        Information about position correspondence.
    residues : list
        For each atom from the enzyme, the residue where it is located.

    """
    aligned = getAligned(srcPDB, tarPDB)
    positions = alignedToPositions(aligned)
    coordinates, residues = getCoordinates(tarPDB, positions)
    return np.array(coordinates), positions, residues


def sphere(params, coordinates):
    """
    The function of the sphere used for least squares fitting.
    """
    x0, y0, z0, r0 = params
    return (coordinates[:, 0] - x0) ** 2 + (coordinates[:, 1] - y0) ** 2 + (coordinates[:, 2] - z0) ** 2 - r0 ** 2


def sphereFitting(coordinates):
    """
    Least squares fitting of the sphere.

    Parameters
    ----------
    coordinates : numpy Array of float64
        The coordinates of the atoms.

    Returns
    -------
    tparams[0] : list
        Parameters of the sphere: x0, y0, z0, r0.

    """
    xm = np.mean(coordinates[:, 0])
    ym = np.mean(coordinates[:, 1])
    zm = np.mean(coordinates[:, 2])
    rm = np.sqrt((coordinates[0, 0] - xm) ** 2 + (coordinates[0, 1] - ym) ** 2 + (coordinates[0, 2] - zm) ** 2)
    
    tparams = leastsq(sphere, [xm, ym, zm, rm], coordinates)
    return tparams[0]


def intersectionCircle(enzymeSphere, activeSphere):
    """
    Calculate the surface area and open angle of the active site

    Parameters
    ----------
    enzymeSphere : list
        The parameters of the enzyme sphere.
    activeSphere : list
        The parameters of the active site sphere.

    Returns
    -------
    sIntersection : float64
        The surface area of the active site on the enzyme.
    alpha : float64
        The opening angles of the active site on the enzyme.

    """
    x0, y0, z0, r0 = enzymeSphere
    x1, y1, z1, r1 = activeSphere
    
    d = np.sqrt((x0 - x1) ** 2 + (y0 - y1) ** 2 + (z0 - z1) ** 2)
    x2 = 2 * (x1 - x0)
    y2 = 2 * (y1 - y0)
    z2 = 2 * (z1 - z0)
    r2Squared = d ** 2 - r1 ** 2 + x1 ** 2 - x0 ** 2 + y1 ** 2 - y0 ** 2 + z1 ** 2 - z0 ** 2
    normal = np.array([x2, y2, z2])
    vector = np.array([x1, y1, z1 - r2Squared / (2 * (z1 - z0))])
    D = np.abs(np.dot(normal, vector)) / np.linalg.norm(normal)
    rIntersection = np.sqrt(np.linalg.norm(vector) - np.linalg.norm(D))
    sIntersection = 2 * np.pi * rIntersection
    alpha = np.arctan(rIntersection / D)
    
    return sIntersection, alpha


def enzyme(srcPDB, tarPDB, searchString):
    """
    Calculate relevant parameters of the active site on a given enzyme

    Parameters
    ----------
    srcPDB : string
        The filename of the PDB file recording the enzyme.
    tarPDB : string
        The filename of the PDB file recording the entire chain (E1-linker-CRY2-CRY2-linker-E2).
    searchString : string
        A string recording the name of residue and id with respect to the enzyme of the active site on the enzyme.

    Returns
    -------
    intersectionCircle(enzymeSphere, activeSphere)
        See the description of intersectionCircle().
    enzymeSphere : list
        The parameters of the enzyme sphere.
    activeSphere : list
        The parameters of the active site sphere.

    """
    coordinates, positions, residues = getAfterCoordinates(srcPDB, tarPDB)
    
    initialResidueID = residues[0].get_id()[1]
    offset = int(searchString[3:]) - initialResidueID
    
    i = 0
    while residues[i].get_id()[1] < initialResidueID + offset:
        i += 1
        
    first = i
    last = -1
    resID = residues[i].get_id()[1]
    
    if i != len(residues):
        last = i + 1
        while last <= len(residues) and residues[last].get_id()[1] == resID:
            last += 1
        last -= 1
    else:
        last = i
    
    global premier, dernier
    premier, dernier = first, last
    enzymeSphere = sphereFitting(coordinates)
    activeSphere = sphereFitting(coordinates[first:last + 1])
    
    return intersectionCircle(enzymeSphere, activeSphere), enzymeSphere, activeSphere


def calculateActiveParam(srcPDB1, srcPDB2, tarPDB, json1, json2):
    """
    Calculate relevant parameters of the active sites on the two enzymes.

    Parameters
    ----------
    srcPDB1 : string
        The filename of the PDB file recording enzyme 1.
    srcPDB2 : string
        The filename of the PDB file recording enzyme 2.
    tarPDB : string
        The filename of the PDB file recording the entire chain (E1-linker-CRY2-CRY2-linker-E2).
    json1 : string
        The filename of the json file recording the name of residue and id with respect to the enzyme of the active site on enzyme 1.
    json2 : string
        The filename of the json file recording the name of residue and id with respect to the enzyme of the active site on enzyme 2.

    Returns
    -------
    sIntersection1 : float64
        The surface area of the active site on enzyme 1.
    alpha1 : float64
        The opening angle of the active site on enzyme 1.
    r1 : float64
        The radius of radius 1.
    sIntersection2 : float64
        The surface area of the active site on enzyme 2.
    alpha2 : float64
        The opening angle of the active site on enzyme 2.
    r2 : float64
        The radius of radius 2.
    d : float64
        The distance between the active sites of enzyme 1 and 2.

    """
    sIntersectionAndAlpha1, enzymeSphere1, activeSphere1= enzyme(srcPDB1, tarPDB, json1)
    sIntersectionAndAlpha2, enzymeSphere2, activeSphere2= enzyme(srcPDB2, tarPDB, json2)
    sIntersection1, alpha1 = sIntersectionAndAlpha1
    sIntersection2, alpha2 = sIntersectionAndAlpha2
    r1 = enzymeSphere1[3]
    r2 = enzymeSphere2[3]
    d = np.sqrt((activeSphere1[0] - activeSphere2[0]) ** 2 + (activeSphere1[1] - activeSphere2[1]) ** 2 + (activeSphere1[2] - activeSphere2[2]) ** 2)
    return sIntersection1, alpha1, r1, sIntersection2, alpha2, r2, d