#!/bin/bash

##################################################
### Printing usage help

print_help()
{
cat << EOF 1>&2

$0 parameters:

  Required:
    -i    path to input file in PDB format
  
  Optional:
    -p    flag to automatically open pymol and load the input PDB file and the produced sript
    -g    residue groups description
    -f    face coloring mode
    -s    selections coloring mode
    -n    output names prefix

  Other:
    -h    show this message and exit
    -e    show more help and exit

EOF
}

print_more_help()
{
cat << EOF 1>&2

What this script does:
  It runs CAD-score voroprot2 program to produce a PyMol (http://pymol.org/) API script
  that draws Voronoi contact faces and, optionally, selects the contacting residues.
  
Running example:
  ./Voroprot2_print_inter_chain_interfaces.bash -i ~/Downloads/2ZSK.pdb -f residue_type -s residue_type -g "(A37-A44)(B37-B44)" -p
  
General notes:
  To run this script you need 'voroprot2' in your binary bath or in the same directory as this script.
  If '-g' option is not provided, inter-chain contacts are considered.
  If '-s' option is not provided, selections are not created.
  To use '-p' option you need 'pymol' in your binary path.

Describing residue groups:
  For example, the string "(A1-A99,B130-B170)(C15)(D)(3-81)" describes
  the contacts between the following four groups of residues:
    Group (A1-A99, B1-B99) stands for the residues from 1 to 99 in chain A and the residues from 130 to 170 in the chain B
    Group (C15, C45) stands for the residues 15 and 45 in the chain C
    Group (D) stands for all the residues in the chain D
    Group (3-81) stands for the residues from 3 to 81 in the unnamed chain 

Available coloring modes for faces and selections:
  blank
  residue_hydrophobicity
  residue_type
  atom_type
  residue_id
  atom_id
ContactAccepterForInterInterval
EOF
}

##################################################
### Reading and checking arguments

SCRIPT_DIRECTORY=$(dirname $0)
VOROPROT_NAME="voroprot2"
VOROPROT="$SCRIPT_DIRECTORY/$VOROPROT_NAME"
if [ ! -f "$VOROPROT" ]
then
  if which $VOROPROT_NAME &> /dev/null
  then
    VOROPROT=$VOROPROT_NAME
  else
    echo "Fatal error: '$VOROPROT_NAME' executable not found" 1>&2
    exit 1
  fi
fi

INPUT_FILE=""
FACE_COLORING_MODE=""
SELECTION_COLORING_MODE=""
COMBINED_RESIDUE_CONTACTS_FILE=""
SPECIFIC_CONTACT_TYPE=""
BINARY_COLORING=""
TRANSPARENT_MAGENTA=""
RESIDUE_GROUPS=""
DRAW_OUTLINE=""
DRAW_INSIDES=""
OUTPUT_NAMES_PREFIX=""
OPEN_IN_PYMOL=false

while getopts "hei:f:s:c:d:vwg:oun:p" OPTION
do
  case $OPTION in
    h)
      print_help
      exit 0
      ;;
    e)
      print_help
      print_more_help
      exit 0
      ;;
    i)
      INPUT_FILE=$OPTARG
      ;;
    f)
      FACE_COLORING_MODE="--face-coloring "$OPTARG
      ;;
    s)
      SELECTION_COLORING_MODE="--selection-coloring "$OPTARG
      ;;
    c)
      COMBINED_RESIDUE_CONTACTS_FILE=$OPTARG
      ;;
    d)
      SPECIFIC_CONTACT_TYPE="--specific-contact-type "$OPTARG
      ;;
    v)
      BINARY_COLORING="--binary-coloring"
      ;;
    w)
      TRANSPARENT_MAGENTA="--transparent-magenta"
      ;;
    g)
      RESIDUE_GROUPS="--groups "$OPTARG
      ;;
    o)
      DRAW_OUTLINE="--outline"
      ;;
    u)
      DRAW_INSIDES="--insides"
      ;;
    n)
      OUTPUT_NAMES_PREFIX="--output-names-prefix "$OPTARG
      ;;
    p)
      OPEN_IN_PYMOL=true
      ;;
    ?)
      exit 1
      ;;
  esac
done

if [ -z "$INPUT_FILE" ]
then
  print_help
  exit 1
fi

if [ ! -f "$INPUT_FILE" ]
then
  echo "Fatal error: input file \"$INPUT_FILE\" does not exist" 1>&2
  exit 1
fi

##################################################
### Generating graphics

TMP_DIR=$(mktemp -d)

SCRIPT_FILE="$TMP_DIR/script.py"

if [ -z "$COMBINED_RESIDUE_CONTACTS_FILE" ]
then
  $VOROPROT --mode collect-atoms < "$INPUT_FILE" | $VOROPROT --mode print-inter-chain-interface-graphics $FACE_COLORING_MODE $SELECTION_COLORING_MODE $RESIDUE_GROUPS $OUTPUT_NAMES_PREFIX $DRAW_OUTLINE $DRAW_INSIDES $SPECIFIC_CONTACT_TYPE> "$SCRIPT_FILE"
else
  ( $VOROPROT --mode collect-atoms < "$INPUT_FILE" ; cat "$COMBINED_RESIDUE_CONTACTS_FILE" ) | $VOROPROT --mode print-inter-chain-interface-graphics --face-coloring inter_residue_contact_scores $SELECTION_COLORING_MODE $RESIDUE_GROUPS $OUTPUT_NAMES_PREFIX $DRAW_OUTLINE $DRAW_INSIDES $SPECIFIC_CONTACT_TYPE $TRANSPARENT_MAGENTA $BINARY_COLORING > "$SCRIPT_FILE"
fi

if [ -s "$SCRIPT_FILE" ]
then
  if $OPEN_IN_PYMOL
  then
    if which pymol &> /dev/null
    then
      pymol "$INPUT_FILE" "$SCRIPT_FILE"
    else
      echo "'pymol' executable not found" 1>&2
    fi
  else
    cat $SCRIPT_FILE
  fi
else
  echo "PyMol script file was not produced" 1>&2
fi

rm -r "$TMP_DIR"
