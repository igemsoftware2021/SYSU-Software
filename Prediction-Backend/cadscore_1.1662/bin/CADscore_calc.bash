#!/bin/bash

set +e

##################################################
### Printing usage help

print_help()
{
cat << EOF 1>&2

$0 parameters:

  Required:
    -D    path to writable database directory (may be not existing)
    -t    path to target file in PDB format
    -m    path to model file in PDB format

  Optional (basic):
    -l    flag to include heteroatoms
    -c    flag to consider only inter-chain contacts
    -z    flag to consider only inter-chain interface zone contacts
    -q    flag to try to rearrange chain names for best possible scores
    -g    flag to use TM-score
    -i    inter-interval contacts specification
    -n    flag to turn on special treatment for nucleic acids

  Optional (advanced):  
    -a    flag to compute atomic global scores
    -r    flag to reset chain names to 'A', 'B', 'C', etc.
    -u    flag to disable model atoms filtering by target atoms
    -b    flag to allow unmatched residue names after filtering
    -s    flag to print summary to standard output
    -y    flag to generate more detailed summary
    -x    flag to delete non-summary data calculated for model
    -j    flag to turn off thread-safe mode
    -v    path to atomic radii files directory
    -e    extra command to produce additional global scores

  Other:
    -h    show this message and exit


A brief tutorial is available from CAD-score site:

  https://bitbucket.org/kliment/cadscore/wiki/Home

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

DATABASE=""
TARGET_FILE=""
MODEL_FILE=""
HETATM_FLAG=""
INTER_CHAIN_FLAG=""
INTERFACE_ZONE_FLAG=""
QUATERNARY_CHAINS_RENAMING=false
USE_TMSCORE=false
INTER_INTERVAL_OPTION=""
USE_ATOMIC_CADSCORE=false
RESETTING_CHAIN_NAMES=""
DISABLE_MODEL_ATOMS_FILTERING=false
ALLOW_UNMATCHED_RESIDUE_NAMES_IN_FILTERING=""
NUCLEIC_ACIDS_MODE=false
GLOBAL_SCORES_CATEGORIES="--categories AA,AM,AS,AW,MA,MM,MS,MW,SA,SM,SS,SW"
PRINT_SUMMARY_TO_STDOUT=false
FULL_GLOBAL_SCORES=false
DELETE_DETAILED_MODEL_DATA=false
THREAD_SAFE_ON=true
RADII_OPTION=""
EXTRA_COMMAND=""

while getopts "hD:t:m:lczqgi:arubnsyxjv:e:" OPTION
do
  case $OPTION in
    h)
      print_help
      exit 0
      ;;
    D)
      DATABASE=$OPTARG
      ;;
    t)
      TARGET_FILE=$OPTARG
      ;;
    m)
      MODEL_FILE=$OPTARG
      ;;
    l)
      HETATM_FLAG="--HETATM"
      ;;
    c)
      INTER_CHAIN_FLAG="--inter-chain"
      ;;
    z)
      INTERFACE_ZONE_FLAG="--interface-zone"
      ;;
    q)
      QUATERNARY_CHAINS_RENAMING=true
      ;;
    g)
      USE_TMSCORE=true
      ;;
    i)
      INTER_INTERVAL_OPTION="--inter-interval "$OPTARG
      ;;
    a)
      USE_ATOMIC_CADSCORE=true
      ;;
    r)
      RESETTING_CHAIN_NAMES="--auto-rename-chains"
      ;;
    u)
      DISABLE_MODEL_ATOMS_FILTERING=true
      ;;
    b)
      ALLOW_UNMATCHED_RESIDUE_NAMES_IN_FILTERING="--allow-unmatched-residue-names"
      ;;
    n)
      NUCLEIC_ACIDS_MODE=true
      GLOBAL_SCORES_CATEGORIES=$GLOBAL_SCORES_CATEGORIES",na_stacking,na_stacking_down,na_stacking_up,na_siding"
      ;;
    s)
      PRINT_SUMMARY_TO_STDOUT=true
      ;;
    y)
      FULL_GLOBAL_SCORES=true
      ;;
    x)
      DELETE_DETAILED_MODEL_DATA=true
      ;;
    j)
      THREAD_SAFE_ON=false
      ;;
    v)
      RADII_OPTION="--radius-classes $OPTARG/vdwr_classes --radius-members $OPTARG/vdwr_members"
      ;;
    e)
      EXTRA_COMMAND=$OPTARG
      ;;
    ?)
      exit 1
      ;;
  esac
done

if [ -z "$DATABASE" ] || [ -z "$TARGET_FILE" ] || [ -z "$MODEL_FILE" ]
then
  print_help
  exit 1
fi

if [ -f "$TARGET_FILE" ]
then
  TARGET_NAME=$(basename $TARGET_FILE)
else
  echo "Fatal error: target file \"$TARGET_FILE\" does not exist" 1>&2
  exit 1
fi

if [ -f "$MODEL_FILE" ]
then
  MODEL_NAME=$(basename $MODEL_FILE)
else
  echo "Fatal error: model file \"$MODEL_FILE\" does not exist" 1>&2
  exit 1
fi

##################################################
### Preparing and checking environment

VERSION_STRING_FILE="$DATABASE/version"
TARGETS_DIR="$DATABASE/targets"
TARGET_DIR="$TARGETS_DIR/$TARGET_NAME"
TARGET_MUTEX_END="$TARGET_DIR/mutex_closed"
TARGET_PARAMETERS_FILE="$TARGET_DIR/parameters"
TARGET_ATOMS_FILE="$TARGET_DIR/atoms"
TARGET_INTER_ATOM_CONTACTS_FILE="$TARGET_DIR/inter_atom_contacts"
TARGET_RESIDUE_IDS_FILE="$TARGET_DIR/residue_ids"
TARGET_INTER_RESIDUE_CONTACTS_FILE="$TARGET_DIR/inter_residue_contacts"
MODELS_DIR="$TARGET_DIR/models"
MODEL_DIR="$MODELS_DIR/$MODEL_NAME"
MODEL_ATOMS_FILE="$MODEL_DIR/atoms"
MODEL_FILTERED_ATOMS_FILE="$MODEL_DIR/filtered_atoms"
MODEL_INTER_ATOM_CONTACTS_FILE="$MODEL_DIR/inter_atom_contacts"
MODEL_RESIDUE_IDS_FILE="$MODEL_DIR/residue_ids"
MODEL_INTER_RESIDUE_CONTACTS_FILE="$MODEL_DIR/inter_residue_contacts"
COMBINED_INTER_RESIDUE_CONTACTS_FILE="$MODEL_DIR/combined_inter_residue_contacts"
CAD_PROFILE_FILE="$MODEL_DIR/cad_profile"
CAD_GLOBAL_SCORES_FILE="$MODEL_DIR/cad_global_scores"
CAD_SIZE_SCORES_FILE="$MODEL_DIR/cad_size_scores"
CAD_ATOMIC_GLOBAL_SCORES_FILE="$MODEL_DIR/cad_atomic_global_scores"
TMSCORE_PROFILE_FILE="$MODEL_DIR/tmscore_profile"
TMSCORE_GLOBAL_SCORES_FILE="$MODEL_DIR/tmscore_global_scores"
EXTRA_COMMAND_GLOBAL_SCORES_FILE="$MODEL_DIR/extra_command_global_scores"
SUMMARY_FILE="$MODEL_DIR/summary"

VERSION_STRING=$($VOROPROT --version | tr -d '\n')

TARGET_PARAMETERS="$HETATM_FLAG $RADII_OPTION $RESETTING_CHAIN_NAMES $INTER_CHAIN_FLAG $INTERFACE_ZONE_FLAG $INTER_INTERVAL_OPTION $NUCLEIC_ACIDS_MODE"

mkdir -p $DATABASE
if [ ! -d "$DATABASE" ] ; then echo "Fatal error: could not create database directory ($DATABASE)" 1>&2 ; exit 1 ; fi

if [ -f "$VERSION_STRING_FILE" ]
then
  CURRENT_VERSION_STRING=$(< $VERSION_STRING_FILE)
  if [ "$VERSION_STRING" != "$CURRENT_VERSION_STRING" ]
  then
    echo "Fatal error: running software version ($VERSION_STRING) is not equal to the version that the database was initialized with ($CURRENT_VERSION_STRING)" 1>&2
    exit 1
  fi
else
  echo -n "$VERSION_STRING" > $VERSION_STRING_FILE
fi

##################################################
### Preprocessing target

mkdir -p $TARGETS_DIR

if mkdir $TARGET_DIR &> /dev/null
then
  echo -n "$TARGET_PARAMETERS" > $TARGET_PARAMETERS_FILE
  
  if [ ! -f $TARGET_ATOMS_FILE ] ; then cat $TARGET_FILE | $VOROPROT --mode collect-atoms $HETATM_FLAG $RADII_OPTION $RESETTING_CHAIN_NAMES > $TARGET_ATOMS_FILE ; fi
  if [ -s "$TARGET_ATOMS_FILE" ] && [ ! -f $TARGET_INTER_ATOM_CONTACTS_FILE ] ; then cat $TARGET_ATOMS_FILE | $VOROPROT --mode calc-inter-atom-contacts > $TARGET_INTER_ATOM_CONTACTS_FILE ; fi
  if [ -s "$TARGET_INTER_ATOM_CONTACTS_FILE" ] && [ ! -f $TARGET_RESIDUE_IDS_FILE ] ; then cat $TARGET_INTER_ATOM_CONTACTS_FILE | $VOROPROT --mode collect-residue-ids  > $TARGET_RESIDUE_IDS_FILE ; fi
  if [ -s "$TARGET_INTER_ATOM_CONTACTS_FILE" ] && [ ! -f $TARGET_INTER_RESIDUE_CONTACTS_FILE ]
  then
    if $NUCLEIC_ACIDS_MODE
    then
      (cat $TARGET_ATOMS_FILE ; cat $TARGET_INTER_ATOM_CONTACTS_FILE | $VOROPROT --mode calc-inter-residue-contacts $INTER_CHAIN_FLAG $INTERFACE_ZONE_FLAG $INTER_INTERVAL_OPTION) | $VOROPROT --mode categorize-inter-nucleotide-side-chain-contacts > $TARGET_INTER_RESIDUE_CONTACTS_FILE
    else
      cat $TARGET_INTER_ATOM_CONTACTS_FILE | $VOROPROT --mode calc-inter-residue-contacts $INTER_CHAIN_FLAG $INTERFACE_ZONE_FLAG $INTER_INTERVAL_OPTION > $TARGET_INTER_RESIDUE_CONTACTS_FILE
    fi
  fi

  true > $TARGET_MUTEX_END
  if [ ! -f "$TARGET_MUTEX_END" ] ; then echo "Fatal error: could not create file ($TARGET_MUTEX_END)" 1>&2 ; exit 1 ; fi

else
  if [ ! -d "$TARGET_DIR" ] ; then echo "Fatal error: could not create target directory ($TARGET_DIR)" 1>&2 ; exit 1 ; fi
  
  if $THREAD_SAFE_ON
  then
    if [ ! -f "$TARGET_MUTEX_END" ]
    then
      TIMEOUT="300"
      WAITING_START_TIME=$(date +%s)
      while [ ! -f "$TARGET_MUTEX_END" ]
      do
  	    CURRENT_TIME=$(date +%s)
  	    WAITING_TIME=$((CURRENT_TIME-WAITING_START_TIME))
  	    if [ "$WAITING_TIME" -gt "$TIMEOUT" ]
  	    then
  	      echo "Fatal error: timeout expired when waiting for target processing" 1>&2
  	      exit 1
  	    fi
      done
    fi
  fi
  
  CURRENT_TARGET_PARAMETERS=$(< $TARGET_PARAMETERS_FILE)
  if [ "$TARGET_PARAMETERS" != "$CURRENT_TARGET_PARAMETERS" ]
  then
    echo "Fatal error: provided script parameters do not match the previous parameters that were used with the same target name and the same database path" 1>&2
    exit 1
  fi
fi

if [ ! -s "$TARGET_ATOMS_FILE" ] ; then echo "Fatal error: no atoms in the target" 1>&2 ; exit 1 ; fi
if [ ! -s "$TARGET_INTER_ATOM_CONTACTS_FILE" ] ; then echo "Fatal error: no inter-atom contacts in the target" 1>&2 ; exit 1 ; fi
if [ ! -s "$TARGET_RESIDUE_IDS_FILE" ] ; then echo "Fatal error: no residues in the target" 1>&2 ; exit 1 ; fi
if [ ! -s "$TARGET_INTER_RESIDUE_CONTACTS_FILE" ] ; then echo "Fatal error: no inter-residue contacts in the target" 1>&2 ; exit 1 ; fi

##################################################
### Preprocessing model

mkdir -p $MODEL_DIR

test -f $MODEL_ATOMS_FILE || cat $MODEL_FILE | $VOROPROT --mode collect-atoms $HETATM_FLAG $RADII_OPTION $RESETTING_CHAIN_NAMES > $MODEL_ATOMS_FILE
if [ ! -s "$MODEL_ATOMS_FILE" ] ; then echo "Fatal error: no atoms in the model" 1>&2 ; exit 1 ; fi

if $DISABLE_MODEL_ATOMS_FILTERING
then
  MODEL_FILTERED_ATOMS_FILE=$MODEL_ATOMS_FILE
fi

test -f $MODEL_FILTERED_ATOMS_FILE || (cat $MODEL_ATOMS_FILE ; cat $TARGET_ATOMS_FILE) | $VOROPROT --mode filter-atoms-by-target $ALLOW_UNMATCHED_RESIDUE_NAMES_IN_FILTERING > $MODEL_FILTERED_ATOMS_FILE
if [ ! -s "$MODEL_FILTERED_ATOMS_FILE" ] ; then echo "Fatal error: no atoms left in the model after filtering by target" 1>&2 ; exit 1 ; fi

test -f $MODEL_INTER_ATOM_CONTACTS_FILE || cat $MODEL_FILTERED_ATOMS_FILE | $VOROPROT --mode calc-inter-atom-contacts > $MODEL_INTER_ATOM_CONTACTS_FILE
if [ ! -s "$MODEL_INTER_ATOM_CONTACTS_FILE" ] ; then echo "Fatal error: no inter-atom contacts in the model" 1>&2 ; exit 1 ; fi

test -f $MODEL_RESIDUE_IDS_FILE || cat $MODEL_INTER_ATOM_CONTACTS_FILE | $VOROPROT --mode collect-residue-ids  > $MODEL_RESIDUE_IDS_FILE
if [ ! -s "$MODEL_RESIDUE_IDS_FILE" ] ; then echo "Fatal error: no filtered residues in the model" 1>&2 ; exit 1 ; fi
	
if [ ! -f $MODEL_INTER_RESIDUE_CONTACTS_FILE ]
then
  if $NUCLEIC_ACIDS_MODE
  then
  	(cat $MODEL_FILTERED_ATOMS_FILE ; cat $MODEL_INTER_ATOM_CONTACTS_FILE | $VOROPROT --mode calc-inter-residue-contacts $INTER_CHAIN_FLAG $INTERFACE_ZONE_FLAG $INTER_INTERVAL_OPTION) | $VOROPROT --mode categorize-inter-nucleotide-side-chain-contacts > $MODEL_INTER_RESIDUE_CONTACTS_FILE
  else
    cat $MODEL_INTER_ATOM_CONTACTS_FILE | $VOROPROT --mode calc-inter-residue-contacts $INTER_CHAIN_FLAG $INTERFACE_ZONE_FLAG $INTER_INTERVAL_OPTION > $MODEL_INTER_RESIDUE_CONTACTS_FILE
  fi
fi
if [ ! -s "$MODEL_INTER_RESIDUE_CONTACTS_FILE" ] ; then echo "Fatal error: no inter-residue contacts in the model" 1>&2 ; exit 1 ; fi

##################################################
### Comparing target and model

if $QUATERNARY_CHAINS_RENAMING
then
  test -f $COMBINED_INTER_RESIDUE_CONTACTS_FILE || cat $TARGET_INTER_RESIDUE_CONTACTS_FILE $MODEL_INTER_RESIDUE_CONTACTS_FILE $TARGET_RESIDUE_IDS_FILE | $VOROPROT --mode calc-combined-inter-residue-contacts --optimally-rename-chains > $COMBINED_INTER_RESIDUE_CONTACTS_FILE
else
  test -f $COMBINED_INTER_RESIDUE_CONTACTS_FILE || cat $TARGET_INTER_RESIDUE_CONTACTS_FILE $MODEL_INTER_RESIDUE_CONTACTS_FILE | $VOROPROT --mode calc-combined-inter-residue-contacts > $COMBINED_INTER_RESIDUE_CONTACTS_FILE
fi

if [ ! -s "$COMBINED_INTER_RESIDUE_CONTACTS_FILE" ] ; then echo "Fatal error: combined inter-residue contacts file is empty" 1>&2 ; exit 1 ; fi
	
test -f $CAD_PROFILE_FILE || cat $COMBINED_INTER_RESIDUE_CONTACTS_FILE $TARGET_RESIDUE_IDS_FILE | $VOROPROT --mode calc-CAD-profile > $CAD_PROFILE_FILE
if [ ! -s "$CAD_PROFILE_FILE" ] ; then echo "Fatal error: CAD profile file is empty" 1>&2 ; exit 1 ; fi

test -f $CAD_GLOBAL_SCORES_FILE || cat $CAD_PROFILE_FILE | $VOROPROT --mode calc-CAD-global-scores $GLOBAL_SCORES_CATEGORIES > $CAD_GLOBAL_SCORES_FILE
if [ ! -s "$CAD_GLOBAL_SCORES_FILE" ] ; then echo "Fatal error: CAD global scores file is empty" 1>&2 ; exit 1 ; fi
	
test -f $CAD_SIZE_SCORES_FILE || cat $CAD_PROFILE_FILE $TARGET_RESIDUE_IDS_FILE $MODEL_RESIDUE_IDS_FILE | $VOROPROT --mode calc-CAD-size-scores > $CAD_SIZE_SCORES_FILE
if [ ! -s "$CAD_SIZE_SCORES_FILE" ] ; then echo "Fatal error: CAD size scores file is empty" 1>&2 ; exit 1 ; fi

if $USE_ATOMIC_CADSCORE
then
  test -f $CAD_ATOMIC_GLOBAL_SCORES_FILE || cat $TARGET_INTER_ATOM_CONTACTS_FILE $MODEL_INTER_ATOM_CONTACTS_FILE | $VOROPROT --mode calc-inter-atom-CAD-score --global $INTER_CHAIN_FLAG > $CAD_ATOMIC_GLOBAL_SCORES_FILE
  if [ ! -s "$CAD_ATOMIC_GLOBAL_SCORES_FILE" ] ; then echo "Fatal error: CAD atomic global scores file is empty" 1>&2 ; exit 1 ; fi
fi

if $USE_TMSCORE
then
  TMSCORE_CALC_NAME="TMscore_calc.bash"
  TMSCORE_CALC="$SCRIPT_DIRECTORY/$TMSCORE_CALC_NAME"
  if [ ! -f "$TMSCORE_CALC" ]
  then
    if which $TMSCORE_CALC_NAME &> /dev/null
    then
      TMSCORE_CALC=$TMSCORE_CALC_NAME
    else
      echo "Fatal error: '$TMSCORE_CALC_NAME' script not found" 1>&2
      exit 1
    fi
  fi
  test -f $TMSCORE_GLOBAL_SCORES_FILE || $TMSCORE_CALC -m $MODEL_FILE -t $TARGET_FILE -p $TMSCORE_PROFILE_FILE -s $TMSCORE_GLOBAL_SCORES_FILE
  if [ ! -s "$TMSCORE_GLOBAL_SCORES_FILE" ] ; then echo "Fatal error: TM-score scores file is empty" 1>&2 ; exit 1 ; fi
fi

if [ -n "$EXTRA_COMMAND" ]
then
  test -f $EXTRA_COMMAND_GLOBAL_SCORES_FILE || $EXTRA_COMMAND $TARGET_FILE $MODEL_FILE > $EXTRA_COMMAND_GLOBAL_SCORES_FILE
  if [ ! -s "$EXTRA_COMMAND_GLOBAL_SCORES_FILE" ] ; then echo "Fatal error: extra command ($EXTRA_COMMAND) scores file is empty" 1>&2 ; exit 1 ; fi
fi

##################################################
### Writing global results

echo target $TARGET_NAME > $SUMMARY_FILE
echo model $MODEL_NAME >> $SUMMARY_FILE

cat $CAD_SIZE_SCORES_FILE >> $SUMMARY_FILE

if $FULL_GLOBAL_SCORES
then
  cat $CAD_GLOBAL_SCORES_FILE >> $SUMMARY_FILE
else
  cat $CAD_GLOBAL_SCORES_FILE | grep -v "_diff" | grep -v "_ref" | grep -v "W" >> $SUMMARY_FILE
fi

if $USE_ATOMIC_CADSCORE ; then cat $CAD_ATOMIC_GLOBAL_SCORES_FILE >> $SUMMARY_FILE ; fi
	
if $USE_TMSCORE ; then cat $TMSCORE_GLOBAL_SCORES_FILE >> $SUMMARY_FILE ; fi
	
if [ -n "$EXTRA_COMMAND" ] ; then cat $EXTRA_COMMAND_GLOBAL_SCORES_FILE >> $SUMMARY_FILE ; fi
	
##################################################
### Optional finalizing

if $PRINT_SUMMARY_TO_STDOUT
then
  cat $SUMMARY_FILE
fi

if $DELETE_DETAILED_MODEL_DATA
then
  rm -f $MODEL_ATOMS_FILE
  rm -f $MODEL_FILTERED_ATOMS_FILE
  rm -f $MODEL_INTER_ATOM_CONTACTS_FILE
  rm -f $MODEL_INTER_RESIDUE_CONTACTS_FILE
fi
