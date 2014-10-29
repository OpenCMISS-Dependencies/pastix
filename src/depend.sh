#!/bin/sh
CC=$1
  shift 1
DEST="$1"
  shift 1
SRC="$1"
  shift 1
OJTDIR="$1"
  shift 1
INCLUDEDIR="$1"
  shift 1
OBJ=`echo $DEST | sed -e "s@\.d@\.o@"`
MODULE=`dirname $OBJ`
MODULE=`dirname $MODULE`
OBJ=`basename $OBJ`

OBJ=`echo $MODULE/$OJTDIR/$OBJ`
TMP=`echo $DEST | sed -e "s@\.d@\.tmp@"`

# echo "$CC $@ < $SRC"
# cpp $(FLAGS) < $(SRC)
# | on garde que les #include 12345 "fichier"
# | on trie et supprime les doublons
# | on vire les dependances non internes
# ! on vire les truc du type MON_OBJET: <built-in> ou MON_OBJET: <command-line>
$CC $@ -DFORCE_NOMPI < $SRC 2>/dev/null| \
    sed -n 's@^\# *[0-9][0-9]* *"\([^"]*\)".*@'$OBJ': \1@p' | \
    sort | uniq | grep -v ": /" | \
    grep -v ": <"> $DEST; [ -s $DEST ] || rm -f $DEST
#sed -e "s@ mpi.h@ @g" -e "s@ pastix_nompi.h@ common/src/nompi.h@g" $DEST > $TMP;
#mv $TMP $DEST;

