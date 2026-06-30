#!/bin/bash
#variable declaration block
po=$PWD
psc=POSCAR
ptc=POTCAR
potloc=/home/p/peckert/VASP_POTCAR
pbe64path=$potloc/PAW_PBE_64
lda64path=$potloc/PAW_LDA_64
pbe54path=$potloc/PAW_PBE_54
lda54path=$potloc/PAW_LDA_54
uspp_ggapath=$potloc/USPP_GGA
uspp_ldapath=$potloc/USPP_LDA
pw91path=$potloc/PAW_PW91
digit="D_sh"
exp="_(sh)"
if [[ $digit =~ "_(sh)" ]]; then
    echo "$digit is a digit"
else
    echo "oops"
fi
#main script

echo "APGU - Automated POTCAR Generation Utility"

echo "Please select the potential you want to use: 1) PBE (PAW) version 64, most recent 2) LDA (PAW) version 64, most recent 3) GGA (USPP), legacy 4) LDA (USPP), legacy 5) PW91 (PAW), legacy 6) PBE (PAW) 54 7) LDA (PAW) 54"
while read potselect
   do
      if [[ $potselect == "1" ]]; then
         potcase=`echo "pbe64"`
         echo "You have selected PBE PAW-based potentials version 64 (most recent as of 15. Nov 2024)."
         break
      elif [[ $potselect == "2" ]]; then
         potcase=`echo "lda64"`
         echo "You have selected LDA PAW-based potentials version 64 (most recent as of 15. Nov 2024)."
         break
      elif [[ $potselect == "3" ]]; then
         potcase=`echo "gga_uspp"`
         echo "You have selected legacy GGA USPP-based potentials."
         break
      elif [[ $potselect == "4" ]]; then
         potcase=`echo "lda_uspp"`
         echo "You have selected legacy LDA USPP-based potentials."
         break
      elif [[ $potselect == "5" ]]; then
         potcase=`echo "pw91"`
         echo "You have selected legacy PW91 (GGA) PAW-based potentials."
         break
     elif [[ $potselect == "6" ]]; then
         potcase=`echo "pbe54"`
         echo "You have selected PBE PAW-based potentials version 54"
         break
      elif [[ $potselect == "7" ]]; then
         potcase=`echo "lda54"`
         echo "You have selected LDA PAW-based potentials version 54"
         break
      else
         echo "Invalid input, run again!"
      fi
   done

function potcar () {
   declare -A pot_dict=( ["no"]="valence shell" ["_h"]="hard" ["_s"]="soft" ["_sv"]="s and p semicore", ["_pv"]="p semicore" ["_d"]="d semicore" ["_2"]="2 f electrons semicore" ["_3"]="3 f electrons semicore")
   local options=()
   local no=1
   local out=""
   for potential in $(find $2 -type d -regextype posix-extended -regex "${2}/${1}(_.*)?"); do #find available potentials
      if [[ ${potential} != *"GW"* ]];then
         out="$out$no) ${potential##*/} " #add potential option to selection list
         options+=($potential)
         no=$((no+1))
      fi
   done
   if [[ $no > 2 ]]; then #if more than two potentials available
      echo "Which potentials do you want to use for ${1}? Potentials without appendix: valcene shell only; _h hard valence shell; _s soft valence shell; _pv p, _sv p and s, _d d semicore states; _n (number) n f electrons semicore states"
      echo $out
      while read typeselect
         do
            if [[ $typeselect < $no ]]; then
               if [[ ${options[$typeselect-1]##*/} =~ _((sv)|[23shd]|(pv)) ]]; then #get potential type selected
                  echo "${pot_dict[${BASH_REMATCH[0]}]} potential for ${1} chosen"
               else
                  echo "Normal valence shell potential for ${1} chosen"
               fi
               cat ${options[$typeselect-1]}/$ptc >> $PWD/$ptc #add potential to POTCAR
            break
            else
               echo "Invalid input, try again!"
            fi
         done
   else
      cat ${options[0]}/$ptc >> $PWD/$ptc
      if [[ ${options[0]##*/} =~ _((sv)|[23shd]|(pv)) ]]; then
         echo "${pot_dict[${BASH_REMATCH[0]}]} potential for ${1} chosen (only potential available)"
      else
         echo "Normal valence shell potential for ${1} chosen (only potential available)"
      fi
   fi
   echo "--------------------------------------------------------------"
}

if [ -f $PWD/$psc ]; then       #if POSCAR exists
   if [ -f $PWD/$ptc ]; then    #delete POTCAR
      rm $PWD/$ptc
   fi
   numatsp=`awk 'NR==6 {print $0}' $PWD/$psc | awk '{print NF}'` #get number of atoms in POSCAR file (6th line)
   for ((atspc=1; atspc<=$numatsp; ++atspc)); do
      atsel=`echo "\$ $atspc" | tr -d ' '`                      #write $number of atom in atsel
      atspid=`awk "NR==6 {print $atsel}" $PWD/$psc | sed 's/[^a-zA-Z]*//g' `            #get name of atom at pos. atsel in POSCAR line 6
      case $potcase in
         pbe64)
         potcar $atspid $pbe64path
         ;;
         lda64)
         potcar $atspid $lda64path
         ;;
         pbe54)
         potcar $atspid $pbe54path
         ;;
         lda54)
         potcar $atspid $lda54path
         ;;
         gga_uspp)
         potcar $atspid $uspp_ggapath
         ;;
         lda_uspp)
         potcar $atspid $uspp_ldapath
         ;;
         pw91)
         potcar $atspid $pw91path
         ;;
      esac
   done
else
   echo "No POSCAR file found in folder: $PWD"
fi
echo "Finished."
