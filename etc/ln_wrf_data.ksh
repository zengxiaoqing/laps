#!/bin/ksh

DATAROOT=/scratch2/portfolios/BMC/public/data/grids/cwb/cwbwrf
LINKROOT='./cwbwrf'

#[ -L $DATAROOT/wrfout* ] && rm $DATAROOT/wrfout*
if [ -L ${LINKROOT}/wrfout* ];then
   print "### REMOVE EXISTING LINKS ###"
   rm -f ${LINKROOT}/wrfout*
fi

# ===== find the lastest dir =====
NDIR=`ls $DATAROOT | tail -1`
if [ -f $DATAROOT/$NDIR/wrfout* ];then
   print "### LASTEST DIR IS : $NDIR ###"
else
   print "### LASTEST DIR IS EMPTY : $NDIR ###"
   NDIR=`ls $DATAROOT | tail -2 | head -1`
   print "### USE THE PREVIOUS ONE : $NDIR ###"
fi

# ===== link wrfout_d02* to wrfout_d01* =====
for FULLNAME in $DATAROOT/$NDIR/wrfout*
do
  WRFD2=${FULLNAME#$DATAROOT/$NDIR/}
  # echo $WRFD2
  WRFD1=${WRFD2:0:9}1${WRFD2:10:20}
  ln -s $FULLNAME $LINKROOT/$WRFD1
done
