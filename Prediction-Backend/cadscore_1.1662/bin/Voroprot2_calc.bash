#!/bin/bash

##################################################
### Printing usage help

print_help()
{
cat << EOF 1>&2

$0 parameters:

  Required:
    -f    path to input file in PDB format

  Optional:
    -l    flag to include heteroatoms
    -a    flag to print atoms
    -q    flag to print quadruples
    -c    flag to print inter-atom contacts
    -r    flag to print inter-residue contacts
    -j    flag to print inter-atom cells faces areas
    -s    flag to print inter-atom contacts summary
    -t    flag to print inter-atom cells faces summary
    -w    flag to check for any inter-chain contact
    -z    cutoff to check for inter-atom clashes
    -v    path to atomic radii files directory
    -u    run zcat on input file
  
  Other:
    -h    show this message and exit

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

MODEL_FILE=""
SPEC_CAT="cat"
HETATM_FLAG=""
RADII_OPTION=""
PRINT_ATOMS=false
QUADRUPLES=false
INTER_ATOM_CONTACTS=false
INTER_RESIDUE_CONTACTS=false
INTER_ATOM_FACES=false
INTER_ATOM_CONTACTS_SUMMARY=false
INTER_ATOM_FACES_SUMMARY=false
CHECK_FOR_ANY_INTER_CHAIN_CONTACT=false
CUTOFF_FOR_CHECKING_CASHES=""

while getopts "hf:laqcrjstwz:v:u" OPTION
do
  case $OPTION in
    h)
      print_help
      exit 0
      ;;
    f)
      MODEL_FILE=$OPTARG
      ;;
    l)
      HETATM_FLAG="--HETATM"
      ;;
    a)
      PRINT_ATOMS=true
      ;;
    q)
      QUADRUPLES=true
      ;;
    c)
      INTER_ATOM_CONTACTS=true
      ;;
    r)
      INTER_RESIDUE_CONTACTS=true
      ;;
    j)
      INTER_ATOM_FACES=true
      ;;
    s)
      INTER_ATOM_CONTACTS_SUMMARY=true
      ;;
    t)
      INTER_ATOM_FACES_SUMMARY=true
      ;;
    w)
      CHECK_FOR_ANY_INTER_CHAIN_CONTACT=true
      ;;
    z)
      CUTOFF_FOR_CHECKING_CASHES="--cutoff $OPTARG"
      ;;
    v)
      RADII_OPTION="--radius-classes $OPTARG/vdwr_classes --radius-members $OPTARG/vdwr_members"
      ;;
    u)
      SPEC_CAT="zcat"
      ;;
    ?)
      exit 1
      ;;
  esac
done

if [ -z "$MODEL_FILE" ]
then
  print_help
  exit 1
fi

if [ ! -f "$MODEL_FILE" ]
then
  echo "File \"$MODEL_FILE\" does not exist" 1>&2
  exit 1
fi

##################################################
### Performing flagged tasks

if $PRINT_ATOMS
then
  $SPEC_CAT $MODEL_FILE | $VOROPROT --mode collect-atoms $HETATM_FLAG $RADII_OPTION
fi

if $QUADRUPLES
then
  $SPEC_CAT $MODEL_FILE | $VOROPROT --mode collect-atoms $HETATM_FLAG $RADII_OPTION | $VOROPROT --mode calc-quadruples
fi

if $INTER_ATOM_CONTACTS
then
  $SPEC_CAT $MODEL_FILE | $VOROPROT --mode collect-atoms $HETATM_FLAG $RADII_OPTION | $VOROPROT --mode calc-inter-atom-contacts
fi

if $INTER_RESIDUE_CONTACTS
then
  $SPEC_CAT $MODEL_FILE | $VOROPROT --mode collect-atoms $HETATM_FLAG $RADII_OPTION | $VOROPROT --mode calc-inter-atom-contacts | $VOROPROT --mode calc-inter-residue-contacts
fi

if $INTER_ATOM_FACES
then
  $SPEC_CAT $MODEL_FILE | $VOROPROT --mode collect-atoms $HETATM_FLAG $RADII_OPTION | $VOROPROT --mode calc-inter-atom-faces
fi

if $INTER_ATOM_CONTACTS_SUMMARY
then
  $SPEC_CAT $MODEL_FILE | $VOROPROT --mode collect-atoms $HETATM_FLAG $RADII_OPTION | $VOROPROT --mode calc-inter-atom-contacts | $VOROPROT --mode summarize-inter-atom-contacts
fi

if $INTER_ATOM_FACES_SUMMARY
then
  $SPEC_CAT $MODEL_FILE | $VOROPROT --mode collect-atoms $HETATM_FLAG $RADII_OPTION | $VOROPROT --mode calc-inter-atom-faces | $VOROPROT --mode summarize-inter-atom-contacts
fi

if $CHECK_FOR_ANY_INTER_CHAIN_CONTACT
then
  $SPEC_CAT $MODEL_FILE | $VOROPROT --mode collect-atoms $HETATM_FLAG $RADII_OPTION | $VOROPROT --mode check-for-any-inter-chain-contact
fi

if [ -n "$CUTOFF_FOR_CHECKING_CASHES" ]
then
  $SPEC_CAT $MODEL_FILE | $VOROPROT --mode collect-atoms $HETATM_FLAG $RADII_OPTION | $VOROPROT --mode check-for-inter-atom-clashes $CUTOFF_FOR_CHECKING_CASHES
fi
