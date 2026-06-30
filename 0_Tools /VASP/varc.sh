#!/bin/bash

echo "Save OUTCAR,CONTCAR,out and vasprun.xml file, e.g. before restarting, also copies CONTCAR to POSCAR after backup."

echo "Checking whether calculation in $PWD has finished..."

fldr=`basename $PWD`

if [[ -e $PWD/OUTCAR ]]; then
   cfin=`grep "Voluntary context" $PWD/OUTCAR | wc -l`
   if [[ $cfin == 1 ]]; then
      yarch=1
   else
      yarch=0
   fi
else
   echo "OUTCAR file does not exist! Calculation probably not started yet. Re-check."
   exit 1
fi

case $yarch in
   0)
   echo "Calculation not yet finished!"
   exit 1
   ;;
   1)
   echo "Backing up files."
   if [[ ! -e $PWD/OUTCAR.s1 ]]; then
      echo "This was the first calculation run - checking for POSCAR_init."
      cp $PWD/OUTCAR $PWD/OUTCAR.s1
      cp $PWD/out_VASP $PWD/out_VASP.s1
      cp $PWD/CONTCAR $PWD/CONTCAR.s1
      cp $PWD/vasprun.xml $PWD/vasprun.xml.s1
      if [[ ! -e $PWD/POSCAR_init ]]; then
         cp $PWD/POSCAR $PWD/POSCAR_init
      else
         echo "POSCAR_init exists already."
      fi
      cp $PWD/CONTCAR $PWD/POSCAR
   else
      echo "This was a subsequent run, will back up files now."
      ocid=`ls -lat OUTCAR.* | head -n 1 | awk '{print $9}' | awk -F '[.]' '{print $2}' | tr -d [:alpha:]` #get number of backup
      otid=`ls -lat out_VASP.* | head -n 1 | awk '{print $9}' | awk -F '[.]' '{print $2}' | tr -d [:alpha:]`
      ccid=`ls -lat CONTCAR.* | head -n 1 | awk '{print $9}' | awk -F '[.]' '{print $2}' | tr -d [:alpha:]`
      vrid=`ls -lat vasprun.xml.* | head -n 1 | awk '{print $9}' | awk -F '[.]' '{print $3}' | tr -d [:alpha:]`
      nocid=$(($ocid + 1))
            notid=$(($otid + 1))
      nccid=$(($ccid + 1))
      nvrid=$(($vrid + 1))
      cp $PWD/OUTCAR $PWD/OUTCAR.s$nocid
      cp $PWD/out_VASP $PWD/out_VASP.s$notid
      cp $PWD/CONTCAR $PWD/CONTCAR.s$nccid
      cp $PWD/vasprun.xml $PWD/vasprun.xml.s$nvrid
      cp $PWD/CONTCAR $PWD/POSCAR
   fi
   ;;
esac

echo "Finished."
