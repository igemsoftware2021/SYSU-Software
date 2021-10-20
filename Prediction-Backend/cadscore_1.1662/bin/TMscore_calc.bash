#!/bin/bash

##################################################
### Printing usage help

print_help()
{
cat << EOF 1>&2

$0 parameters:

  Required:
    -t    path to target file in PDB format
    -m    path to model file in PDB format
    -p    profile output file
    -s    summary output file

  Other:
    -h    show this message and exit
  
Note: $0 needs TMscore application
      (http://zhanglab.ccmb.med.umich.edu/TM-score/)
      to be available either in the same directory as this script
      or in your system binary path

EOF
}

##################################################
### Reading and checking arguments

SCRIPT_DIRECTORY=$(dirname $0)
TMSCORE_BIN_NAME="TMscore"
TMSCORE_BIN="$SCRIPT_DIRECTORY/$TMSCORE_BIN_NAME"
if [ ! -f "$TMSCORE_BIN" ]
then
  if which $TMSCORE_BIN_NAME &> /dev/null
  then
    TMSCORE_BIN=$TMSCORE_BIN_NAME
  else
    echo "Fatal error: '$TMSCORE_BIN_NAME' executable not found" 1>&2
    exit 1
  fi
fi

TARGET_FILE=""
MODEL_FILE=""
TMSCORE_PROFILE_FILE=""
TMSCORE_SUMMARY_FILE=""

while getopts "ht:m:p:s:" OPTION
do
  case $OPTION in
    h)
      print_help
      exit 0
      ;;
    t)
      TARGET_FILE=$OPTARG
      ;;
    m)
      MODEL_FILE=$OPTARG
      ;;
    p)
      TMSCORE_PROFILE_FILE=$OPTARG
      ;;
    s)
      TMSCORE_SUMMARY_FILE=$OPTARG
      ;;
    ?)
      exit 1
      ;;
  esac
done

if [ -z "$TARGET_FILE" ] || [ -z "$MODEL_FILE" ] || [ -z "$TMSCORE_PROFILE_FILE" ] || [ -z "$TMSCORE_SUMMARY_FILE" ]
then
  print_help
  exit 1
fi

if [ ! -f "$TARGET_FILE" ]
then
  echo "Target file \"$TARGET_FILE\" does not exist" 1>&2
  exit 1
fi

if [ ! -f "$MODEL_FILE" ]
then
  echo "Model file \"$MODEL_FILE\" does not exist" 1>&2
  exit 1
fi

##################################################
### Getting values from TM-score

TEMP_SHORT_TARGET_FILE=$(mktemp)
TEMP_SHORT_MODEL_FILE=$(mktemp)

cp $TARGET_FILE $TEMP_SHORT_TARGET_FILE
cp $MODEL_FILE $TEMP_SHORT_MODEL_FILE

$TMSCORE_BIN $TEMP_SHORT_MODEL_FILE $TEMP_SHORT_TARGET_FILE > $TMSCORE_PROFILE_FILE
rm $TEMP_SHORT_MODEL_FILE $TEMP_SHORT_TARGET_FILE

TM_SCORE=$(cat $TMSCORE_PROFILE_FILE | egrep "TM-score\s*=.*d0" | sed 's/TM-score\s*=\s*\(.*\)\s*(.*/\1/g')
if [ -z "$TM_SCORE" ] ; then TM_SCORE=0 ; fi

TM_SCORE_GDT_TS=$(cat $TMSCORE_PROFILE_FILE | egrep "GDT-TS-score" | sed 's/GDT-TS-score\s*=\s*\(.*\)\s*%(d<1).*/\1/g')
if [ -z "$TM_SCORE_GDT_TS" ] ; then TM_SCORE_GDT_TS=0 ; fi

TM_SCORE_GDT_HA=$(cat $TMSCORE_PROFILE_FILE | egrep "GDT-HA-score" | sed 's/GDT-HA-score\s*=\s*\(.*\)\s*%(d<0\.5).*/\1/g')
if [ -z "$TM_SCORE_GDT_HA" ] ; then TM_SCORE_GDT_HA=0 ; fi

echo "TM_score $TM_SCORE" > $TMSCORE_SUMMARY_FILE
echo "TM_score_GDT_TS $TM_SCORE_GDT_TS" >> $TMSCORE_SUMMARY_FILE
echo "TM_score_GDT_HA $TM_SCORE_GDT_HA" >> $TMSCORE_SUMMARY_FILE
