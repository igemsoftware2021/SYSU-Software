# Overview

## About this project

CAD-score is a command-line software that constructs the Voronoi diagram of protein atom spheres, calculates inter-atom contact surfaces and uses them to compare protein structures of the same sequence.

Please note that this document describes the old, original implementation of CAD-score.

## Differences between the old and the new implementations

The new CAD-score implementation is the "voronota-cadscore"
script in the Voronota package that is available from [https://github.com/kliment-olechnovic/voronota/releases](https://github.com/kliment-olechnovic/voronota/releases).

In the older implementation, contacts are constructed for every atom by subdividing the expanded atom sphere according to the Voronoi neighbors.
In the more recent implementation, contacts are derived directly from the Voronoi faces.
The scores produced by implementations differ, but correlate very strongly.
Note that the new implementation can run in the classic mode with the "--old-regime" flag.

## Links

[Official CAD-score method page.](http://bioinformatics.lt/software/cad-score)

[CAD-score web-interface.](http://bioinformatics.ibt.lt/cad-score)

# Using CAD-score

## Getting the software

To use this standalone CAD-score software, download the latest package
from [https://github.com/kliment-olechnovic/old_cadscore/releases](https://github.com/kliment-olechnovic/old_cadscore/releases).

CAD-score software project consists of:

* "voroprot2" C++ application
* Several Bash scripts that run "voroprot2" to perform various tasks 

Every script in this project can print a list of expected parameters when executed without any arguments or with "-h" flag.

## Preparing the software

A compiled static executable for GNU/Linux is included in the package. However, it is recommended to rebuild the executable using the following command:

    g++ -O3 -o bin/voroprot2 src/*.cpp

Users of very old C++ old compilers may need to provide an additional option:

    g++ -DFOR_OLDER_COMPILERS -O3 -o bin/voroprot2 src/*.cpp

## Basic command-line usage example

Assume that we want to score two protein structure models "[model1.pdb](https://raw.githubusercontent.com/kliment-olechnovic/old_cadscore/master/tests/basic/input/model1)"
and "[model2.pdb](https://raw.githubusercontent.com/kliment-olechnovic/old_cadscore/master/tests/basic/input/model2)"
against the reference structure "[target.pdb](https://raw.githubusercontent.com/kliment-olechnovic/old_cadscore/master/tests/basic/input/target)"
(note that residue sequence, residue numbering and chains naming in the models should be consistent with the target).
Scoring is done by running the following commands:

    CADscore_calc.bash -D /path/to/database -t /path/to/target.pdb -m /path/to/model1.pdb
    CADscore_calc.bash -D /path/to/database -t /path/to/target.pdb -m /path/to/model2.pdb

The comparison data is now contained in the "database" directory. The global scores table is collected from the "database" directory using the following command:

    CADscore_read_global_scores.bash -D /path/to/database

Local contact differences for the "target.pdb" model "model1.pdb" are collected using the "CADscore_read_local_scores.bash":

    CADscore_read_local_scores.bash -D /path/to/database -t target.pdb -m model1.pdb -c AS -w 3

Here "-c AS" means that we are interested in "A-S" contacts and "-w 3" means that that we want each value to be smoothed by window of (3+1+3) positions.

If [TMscore](http://zhanglab.ccmb.med.umich.edu/TM-score/) program is available in your system binary path, you can use "-g" flag to tell "CADscore_calc.bash" to additionally compute TM-score, GDT-TS and GDT-HA global scores.

Note that CAD-score uses file base-names, not file full-paths as identifiers. So, for a single database directory, base-names of target files should be unique. And, for each target in a database, base-names of model files should be unique.

## Evaluation modes

CAD-score evaluation modes are:

* All contacts (default mode): all contacts within the target structure (single-domain, multidomain or multichain)
* Inter-chain contacts: contacts between different chains
* Custom contacts: contacts between user-defined residue ranges 

"Inter-chain" mode is turned on with "-c" flag:

    CADscore_calc.bash -D database_for_inter_chain -t target.pdb -m model1.pdb -c

"Custom contacts" mode is turned on with "-i" option:

    CADscore_calc.bash -D database_for_custom -t target.pdb -m model1.pdb -i "(A1-A99)(A100-A150)"

Custom contacts are described as a list of residue groups. For example, the string "(A1-A99, B130-B170) (C15) (D) (3-81)" forces CAD-score to consider only the contacts between the following four groups of residues:

1. Group (A1-A99, B130-B170) stands for the residues from 1 to 99 in chain A and the residues from 130 to 170 in the chain B
2. Group (C15, C45) stands for the residues 15 and 45 in the chain C
3. Group (D) stands for all the residues in the chain D
4. Group (3-81) stands for the residues from 3 to 81 in the unnamed chain 

## Using CAD-score for clustering

Assume that we want cluster protein structure models contained in a directory. CAD-score can calculate similarity matrices for them and then call R to draw heatmaps and dendrograms:

    CADscore_create_scores_matrices.bash -I models_dir -O output_dir -d

## Note on comparing homo-oligomers

When comparing homo-oligomers, it is not always obvious which chain in the model corresponds to which chain in the target. "-q" option can be used to rearrange the chain names in the homo-oligomer model to get the highest possible CAD-score values:

    CADscore_calc.bash -D database -t target.pdb -m model.pdb -q

## Note on comparing RNA structures

Use "-n" option for distinguishing stacking and pairing interactions in RNA structures.

