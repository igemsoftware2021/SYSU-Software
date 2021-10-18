#!/bin/bash

##################################################
### Printing usage help

print_help()
{
cat << EOF 1>&2

$0 parameters:

  Required:
    -D    path to database directory
    -t    target name in the database
    -m    model name in the database
    -c    contacts category
    -a    flag for non-normalized output

  Optional:
    -w    bluring window size

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

DATABASE=""
TARGET_NAME=""
MODEL_NAME=""
CATEGORY=""
ABSOLUTE_FLAG=""
WINDOW="0"

while getopts "hD:t:m:c:aw:" OPTION
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
      TARGET_NAME=$(basename "$OPTARG")
      ;;
    m)
      MODEL_NAME=$(basename "$OPTARG")
      ;;
    c)
      CATEGORY=$OPTARG
      ;;
    a)
      ABSOLUTE_FLAG="--absolute"
      ;;
    w)
      WINDOW=$OPTARG
      ;;
    ?)
      exit 1
      ;;
  esac
done

if [ -z "$DATABASE" ] || [ -z "$TARGET_NAME" ] || [ -z "$MODEL_NAME" ] || [ -z "$CATEGORY" ]
then
  print_help
  exit 1
fi

##################################################
### Processing CAD profile

CAD_PROFILE_FILE="$DATABASE/targets/$TARGET_NAME/models/$MODEL_NAME/cad_profile"

if [ ! -s "$CAD_PROFILE_FILE" ]
then
  echo "Contacts area difference profile file for target \"$TARGET_NAME\" and model \"$MODEL_NAME\" is not in the database" 1>&2
  exit 1
fi

cat $CAD_PROFILE_FILE | $VOROPROT --mode calc-CAD-local-scores --category $CATEGORY --window $WINDOW $ABSOLUTE_FLAG
