#!/bin/bash
PDIR=../../../out
CUR=`pwd`
DIR=`basename $CUR`
 
WDIR=$PDIR/$DIR
DDIR=Current_Results
cp $DDIR/CHRISTIFIRE*devc.csv $WDIR
cp $DDIR/*git.txt $WDIR

