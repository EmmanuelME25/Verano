#!/bin/bash
##
#
# This programs generates the kpoints in the IBZ
# for BCC (Face Centered Cubic) and the k-Fermi-Surface
# for the Real FS or the ideal Sommerfeld (spherical) FS
# using a refinement scheme
#
#  !* input files
#  !* fort.1
#  !* fort.4
#  !* fort.21
#  !* fort.22
#  !* fort.35
#  !* fort.244
#  !* fort.245
#  !* fort.750
#  !* fort.751
#  !* fort.246
#  !* fort.247
#  !* fort.248
#
##
red='\e[0;31m'
RED='\e[1;31m'
blue='\e[0;34m'
BLUE='\e[1;34m'
cyan='\e[0;36m'
CYAN='\e[1;36m'
GREEN='\e[0;32m'
GRE='\e[1;32m'
YELLOW='\e[1;33m'
MAG='\e[0;35m'
NC='\e[0m' # No Color
##

## debug
function despulga {
Line
printf "\taqui\n"
Line
exit 1
}
## thanks
function gracias {
    printf "\tThanks for using ${cyan}TINIBA${NC}: ${RED}NO WARRANTIES WHATSOEVER\n${NC}"
}
##
function Line {
    printf "\t${cyan}--------------------${NC}\n"
}
##
case=`echo $PWD | awk -F / '{print$NF}'`
exec="$TINIBA/utils/brillouin-zone/ibz/bcc/rbands-fs-bz-bcc-grid-integrate"
script=bands-fs-bz-bcc-grid-integrate.sh
sym=symmetries
## input at will
########### this is really set in responses.sh
if [ "$#" -eq 0 ]
then   # Script needs at least one command-line argument.
    clear
    Line
    printf "\tGenerates IBZ-BCC refined Fermi Surface k-points starting from a Coarse-Grained Fermi Surface\n"
    printf "\tvia bands-fs-bz-bcc-grid-integrat.sh${NC}\n\n"
    sh linea.sh
    flash.sh "Warning"
    printf "\tThe Coarse-Grained Fermi Surface must be previously calculated via\n"
    printf "\t${blue}new-fermi-surface-via-kf.sh${NC}\n\n"
    sh linea.sh
    printf "\tUsage:\n"
    printf "\t$script -N ${red}i${NC} -M ${red}i${NC} -t ${red}r${NC} -b ${red}i${NC} -c ${red}r${NC} -p ${red}r${NC} -e ${red}r${NC} -k ${red}r${NC}\n"
    printf "\n\t where\n"
    printf "\t-N ${red}Fine-Grained Grid${NC}\n"
    printf "\t-M ${red}Coarse-Grained Grid${NC}\n"
    printf "\t-t ${red}radius around CG-FS ${NC}\n"
    printf "\t-b ${red}band${NC}\n"
    printf "\t-c ${red}kc${NC}\n"
    printf "\t-p ${red}percentage${NC}\n"
    printf "\t-e ${red}E_Fermi${NC}\n"
    printf "\t-k ${red}k-centroid for band 36 or 37, otherwise use 0${NC}\n\n"
    printf "\ti->integer, r->real\n"
    Line
    exit 1
fi
# gets options
vnlkss=false
while getopts ":N:M:b:a:c:p:e:k:t:" OPTION
do
    case $OPTION in
        N)
	    N=$OPTARG
	    ;;
        M)
	    Ncg=$OPTARG
	    ;;
        b)
	    band=$OPTARG
	    ;;
        c)
	    kc=$OPTARG
	    ;;
        p)
	    p=$OPTARG
	    ;;
        e)
	    efermi=$OPTARG
	    ;;
        k)
	    kcBand=$OPTARG
	    ;;
        t)
	    tol=$OPTARG
	    ;;
        ?)
	printf "\t${RED}error${NC}\n"
	exit
	;;
    esac
done
###Let the Script Begin
#info
echo $Ncg $N > .info
#exit
como_corri="$script -N $N -M $Ncg -t $tol -b $band -c $kc -p $p -e $efermi -k $kcBand >> como-corri-$script"
date >> como-corri-$script
echo $script -N $N -M $Ncg -t $tol -b $band -c $kc -p $p -e $efermi -k $kcBand >> como-corri-$script
#
sh linea.sh
sh begin.sh
printf "\t$como_corri\n"
sh linea.sh
#getting info
sym=symmetries
#gets the value of a,b,c
sh linea.sh
sh acell.sh
a=`awk '{print $1}' fort.1`;a=`echo "scale=10;$a/1" | bc`;a=`echo "scale=10;$a/1" | bc`
b=`awk '{print $2}' fort.1`;b=`echo "scale=10;$b/1" | bc`;b=`echo "scale=10;$b/1" | bc`
c=`awk '{print $3}' fort.1`;c=`echo "scale=10;$c/1" | bc`;c=`echo "scale=10;$c/1" | bc`
echo $a $b $c > fort.1
eta=`echo "scale=4;$c/$a" | bc`
Line
printf "\tUnit cell: a=$a b=$b c=$c Bohrs\n"
# patito
#set view 90,225,2.5
#Full BZ
echo "a=$a" > $sym/bz-bcc.g
cat $TINIBA/utils/brillouin-zone/ibz/bcc/only-bz-bcc.g >> $sym/bz-bcc.g
#IBZ
echo "a=$a" > $sym/ibz-bcc.g
cat $TINIBA/utils/brillouin-zone/ibz/bcc/only-ibz-bcc.g >> $sym/ibz-bcc.g
#spin
while read value
do
    spin=$value
done < .spininfo
#
grep ecut setUpAbinit_$case.in > .zaz
ecut=`awk '{print $2}' .zaz`
#
if [[ "$band" == "34" || "$band" == "35" ]]
then
    #
    #Extracting nm for grad=-1
    #&          np for grad=+1
    #For bands 34 and 35 nm=np
    cgf=.coarse-grained_$Ncg
    nm=`awk 'FNR==2{print $2}' $cgf`
    np=`awk 'FNR==2{print $3}' $cgf`
    #
    sh linea.sh
    printf "\tFrom $cgf: $Ncg nm=$nm np=$np\n"
    #For bands 34 and 35 nm=np
    finm=$sym/$case.fs-ibz-in-cg-$Ncg-r-$Ncg-k-$nm
    printf "\twe use: $Ncg Coarse-Grained FS-k-points with $nm for grad=-1@\n"
    printf "\t$finm -> fort.244\n"
    cp $finm fort.244
fi
if [[ "$band" == "36" || "$band" == "37" ]]
then
    #
    #Extracting nm for grad=-1
    #&          np for grad=+1
    #For bands 36 and 37 nm!=np
    cgf=.coarse-grained_$Ncg
    total=`awk 'FNR==2{print $2}' $cgf`
    nm=`awk 'FNR==2{print $3}' $cgf`
    np=`awk 'FNR==2{print $4}' $cgf`
    #
    finm=$sym/$case.fs-ibz-in-cg-$Ncg-r-$Ncg-k-$total-1-$nm
    finp=$sym/$case.fs-ibz-in-cg-$Ncg-r-$Ncg-k-$total+1-$np
    printf "\twe use: $Ncg  divisions for $total Coarse-Grained FS-k-points with $nm for grad=-1 & $np for grad=+1 @\n"
    printf "\t$finm -> fort.244\n"
    printf "\t$finp -> fort.245\n"
    cp $finm fort.244
    cp $finp fort.245
fi
#
sh linea.sh
#
sh linea.sh
printf "\t@$script: ${red}For this run the value of Nk=$Ncg is irrelevant${NC}\n"
sh linea.sh
printf "\techo $N $Ncg $a $kc $p $efermi $kcBand $band $nm $np $tol 1 | $exec\n"
sh linea.sh
#     1  2    3  4   5  6       7       8     9   10  11  12
echo $N $Ncg $a $kc $p $efermi $kcBand $band $nm $np $tol 1 | $exec
sh linea.sh
printf "\t$exec ${red}is done${NC}!\n"
sh linea.sh
#
#################Exit-Down
if [ "1" == "2" ]
   then
       sh linea.sh
       flash.sh "@$script:Exit#1\n"
       sh linea.sh
       exit 1
fi
#################Exit-Up
#sort fort.21 & fort.751: valid for all bands--Down
#
sort -u fort.21 > sort.21
nkm=`wc sort.21 | awk '{print $1}'` #N k-points for grad=-1
f1=$sym/$case.kcartesian-grad-1-N-$N-Nk-$nkm
cp sort.21 $f1
sort -u fort.751 > sort.751
nkm=`wc sort.751 | awk '{print $1}'` #N k-points for grad=-1
f4=$sym/$case.selected-grad-1-N-$N-Nk-$nkm
cp sort.751 $f4
cat sort.751 > .auxi
#
#sort fort.21 & fort.751--UP
#
#
#sort fort.22 & fort.750: valid only for band 36 & 37--Down
#
if [[ "$band" == "36" ||  "$band" == "37"  ]]
then
    sort -u fort.22 > sort.22
    nkp=`wc sort.22 | awk '{print $1}'` #N k-points for grad=+1
    f2=$sym/$case.kcartesian-grad+1-N-$N-Nk-$nkp
    cp sort.22 $f2
    sort -u fort.750 > sort.750
    nkp=`wc sort.750 | awk '{print $1}'` #N k-points for grad=+1
    f3=$sym/$case.selected-grad+1-N-$N-Nk-$nkp
    cp sort.750 $f3
    #data for next run of run_tiniba.sh
    cat sort.750 sort.751 > .auxi
fi
#
#sort fort.22 & fort.750--Up
#
nn=`wc .auxi | awk '{print $1}'` #refined k-points for grad=+/- 1
nnk=$case.rkcartesian-$N-$nn
cp .auxi $nnk
fsn=$sym/$case.kcartesian_$nn
mv .auxi $fsn
printf "\t@$script: ${red}$fsn${NC} for refined FS\n"
sh linea.sh
printf "\t@$script: Converting Cartesian ${red}$nnk${NC} to direct lattice kpoints\n"
#
awk '{print $1,$2,$3}' $nnk > .auxfs
#Converts cartesian to direct k-points
#kcart2dlatt.sh -f k-file -o (1,2)[Cartesian->DL,DL->Cartesian]
kcart2dlatt.sh -f .auxfs -o 1
Nkfs=`wc .auxfs-direct | awk '{print $1}'`
#echo $Nkfs > .knfs_$nk
out21d=$case.klist_$Nkfs
mv .auxfs-direct $out21d
#
printf "\tOutput:\n"
if [[ "$band" == "34" || "$band" == "35" ]]
then
    printf "\t       $finm => original $Ncg divisions grad=-1 seed\n"
    printf "\t       $f1 => all cartesian\n"
    printf "\t       $f4 => selected cartesian\n"
    printf "\t       ${red}$nnk${NC} => ${CYAN}refined cartesian used in next run${NC}\n"
fi
if [[ "$band" == "36" || "$band" == "37" ]]
then
    printf "\t       $finm => original $Ncg divisions grad=-1 seed\n"
    printf "\t       $finp => original $Ncg divisions grad=+1 seed\n"
    printf "\t       $f1 => all cartesian\n"
    printf "\t       $f2 => all cartesian\n"
    printf "\t       $f3 => selected cartesian\n"
    printf "\t       $f4 => selected cartesian\n"
    printf "\t       ${red}$nnk${NC} => ${CYAN}refined cartesian used in next run${NC}\n"
fi
################################################1################################################
#gunplot-Down
outdir=plots
if [ ! -d $outdir ]
   then
       mkdir $outdir
fi
goutf=$outdir/band-$band-CG-$Ncg-N-$N-fs.g
printf "\tgnuplot:\n"
printf "\t       gnuplot> load 'symmetries/ibz-bcc.g'\n"
printf "\t       gnuplot> load '$goutf'\n"
#echo "set key at 0,g/2,-g/20" >> $goutf
echo "unset key" > $goutf
echo "set key" >> $goutf
echo "set ticslevel 0" >> $goutf
echo "set view 78,274,1" >> $goutf
echo "set xlabel 'x'" >> $goutf
echo "set ylabel 'y'" >> $goutf
echo "set zlabel 'z'" >> $goutf
sp=.6
cl=.8
echo "sp '$f1' u 1:2:3 w p pt 7 ps $sp lc $cl t 'all:-1'" >> $goutf
if [[ "$band" == "36" ||  "$band" == "37"  ]]
then
    cl=2
    echo "rep '$f2' u 1:2:3 w p pt 7 ps $sp lc $cl t 'all:+1'" >> $goutf
    sp=.8
    cl='rgb "red"'
    echo "rep '$f3' u 1:2:3 w p pt 7 ps $sp lc $cl t 'refined:+1'" >> $goutf
fi
cl='rgb "blue"'
echo "rep '$f4' u 1:2:3 w p pt 7 ps $sp lc $cl t 'refined:-1'" >> $goutf
sp=1
cl='rgb "green"'
echo "rep '$finm' u 1:2:3 w p pt 7 ps $sp lc $cl t 'CG:-1'" >> $goutf
if [[ "$band" == "36" ||  "$band" == "37"  ]]
then
    cl='rgb "brown"'
    echo "rep '$finp' u 1:2:3 w p pt 7 ps $sp lc $cl t 'CG:+1'" >> $goutf
fi
#gunplot-Up
#################Exit-Down
####################################
if [ "1" == "2" ]                  #
   then                            #
       sh linea.sh                 #
       flash.sh "@$script:Exit#2"  #
       sh linea.sh                 #
       exit 1                      #
fi                                 #
####################################
#################Exit-Up
###########################################################################
# we must run run_tiniba.sh to calculate the energy around the refined FS #
###########################################################################
sh linea.sh
printf "\t${CYAN}$out21d => used by run_tiniba.sh for refined $Nkfs k-points in search of FS${NC}\n"
sh linea.sh
#################Exit-Down
####################################
if [ "1" == "2" ]                  #
   then                            #
       sh linea.sh                 #
       flash.sh "@$script:Exit#3"  #
       sh linea.sh                 #
       exit 1                      #
fi                                 #
####################################
#################Exit-Up
######################################################
################################# run_tiniba.sh-Down #
cuale=eigen_$Nkfs\_$ecut-$spin
printf "\tCalculating $cuale\n"
sh linea.sh
if [ ! -f $cuale ]
then
    #### Commands
    nodo=`hostname`
    if [[ "$nodo" == "fat1" || "$nodo" == "fat2" || "$nodo" == "fat3" ]]
    then
	cores=64
    else #hexas
	cores=12
    fi
    if [[ "$nodo" == "viginti01" || "$nodo" == "viginti02" ]]
    then
	cores=40
    fi
    #Number of plane waves
    pw=`grep mpw $case'_check'/$case.out | awk '{print $12}'`
    sh linea.sh
    printf "\t@$script: using $cores cores for SCF\n"
    sh linea.sh
    #here we don't need -p, since these are not the refined FS-k-points
    exectiniba="run_tiniba.sh -r run -k $Nkfs  -N 0 -x 2 -C $cores -P $pw -w -e"
    printf "\t${blue}####### executing for energies only ######${NC}\n"
    printf "\t${red}$exectiniba\n"
    printf "\t${blue}##################################################################${NC}\n"
    # executable
    $exectiniba
    #
    sh linea.sh
    printf "\t${cyan}run_tiniba.sh done!${NC} ${red}back @$script${NC}\n"
    sh linea.sh
else
    flash.sh "case $cuale done!"
    sh linea.sh
fi
################################# run_tiniba.sh-Up ##
#################Exit-Down
####################################
if [ "1" == "2" ]                  #
   then                            #
       sh linea.sh                 #
       flash.sh "@$script:Exit#4"  #
       sh linea.sh                 #
       exit 1                      #
fi                                 #
####################################
#################Exit-Up
#####################################################
########
flash.sh "Selecting next refinement of the FS using $Nkfs k-points"
sh linea.sh
########
########################################################################################################
############################################# selecting the band & pasting the corresponding k-points ##
########################################################################################################
eband=`echo "scale=10;$band+1" | bc`
energy=eigen_$Nkfs\_$ecut-$spin
printf "\t@$script: Using band $band of $energy (=>column $eband) \n"
awk '{print $eband}' $energy > .eaux
echo "awk '{print \$$eband}' $energy > .energy-$band" > .run
chmod +x .run
./.run
out=$sym/$case.kband_$band\_$Nkfs
flash.sh "WARNING"
printf "\tPasting $sym/$case.kcartesian_$Nkfs .energy-$band > $out\n"
sh linea.sh
paste $sym/$case.kcartesian_$Nkfs .energy-$band > $out
outa=$sym/$case.kcartesian_$Nkfs
rm .run .eaux .zaz
sh linea.sh
printf "\t@${red}$out${NC} as input for ${blue}next refined FS${NC} ${cyan}@fort.1${NC}\n"
sh linea.sh
#################Exit-Down
####################################
if [ "1" == "2" ]                  #
   then                            #
       sh linea.sh                 #
       flash.sh "@$script:Exit#5"  #
       sh linea.sh                 #
       exit 1                      #
fi                                 #
####################################
#################Exit-Up
###########################################################################
#we must run bands-fs-bcc-integrate.f90 to generate the new FS k-points #
###########################################################################
#out-directories
sym=symmetries
fs=fermi-surface
if [ ! -d $sym ]
then
    mkdir $sym
fi
if [ ! -d $fs ]
then
    mkdir $fs
fi
sh linea.sh
#
cp $out fort.1 #read by fs-bz-bcc-grid-integrate.f90 with energy(k)
printf "\t@$script: Using $out for band=$band to run:\n"
printf "\techo $N $Nkfs $a $kc $p $efermi $kcBand $band $nm $np $tol 2 | $exec\n"
sh linea.sh
######below command is where the next set of refined k-points are generated#####
echo $N $Nkfs $a $kc $p $efermi $kcBand $band $nm $np $tol 2 | $exec           #
#that runs: data-grad-1+1.sh                                                   #
################################################################################
#
#INFO for PARALELIZE the do loops that generate the refined grid from the previous grid
#@bands-fs-bz-bcc-grid-integrate.f90
#Look for   !!!! @@@@@@@@@@@@PARALELIZE@@@@@@@@@@@
dirpa=paralelize
if [ ! -d $dirpa ]
   then
       mkdir $dirpa
fi
cp fort.21  $dirpa/fort.21-Ncg-$Ncg-N-$N
cp fort.244 $dirpa/fort.244-Ncg-$Ncg-N-$N
sh linea.sh
printf "\tFor paralelize\n"
printf "\t$dirpa/fort.21-Ncg-$Ncg-N-$N\n"
printf "\t$dirpa/fort.244-Ncg-$Ncg-N-$N\n"
sh linea.sh
#
#################Exit-Down
if [ "1" == "2" ]
   then
       sh linea.sh
       flash.sh "@$script:Exit#6"
       #printf "\t${red}@$script: $exec is NOT working${NC}\n"
       #printf "\t${red}@$script: fort.244 disappears ${NC}\n"
       sh linea.sh
       exit 1
fi
#################Exit-Up
#################
#Output
nkfe=`wc fort.4 | awk '{print $1}'` #FULL SET OF k-POINTS such that kmin < k < kmax for the calculation of the FS
out4=$sym/$case.fs-ibz\_$nkfe
#sort fort.4
sort -u fort.4 > sort.4
cp sort.4 $out4
flash.sh "all k-points for FS"
printf "\tOUTPUT: $out4\n"
#
if (( $(echo "$kcBand == .0" | bc -l) ))
then
    kfs=`wc fort.244 | awk '{print $1}'` #number of FS k-points
    ffs=$sym/$case.fs-ibz-in-$kfs
    cp fort.244 $ffs
    kfso=`wc fort.245 | awk '{print $1}'` #number of k-points just above FS
    ffso=$sym/$case.fs-ibz-out-$kfso
    cp fort.245 $ffso
    printf "\tFS->$ffs | Just above FS->$ffso \n"
    sh linea.sh
fi
#
outdir=octant
if [ ! -d $outdir ]
   then
       mkdir $outdir
fi
### @@@@-D
# area-file
af=area-vs-k-fs.dat
if [ ! -f $af ]
then
    touch $af
    echo "#     1N    2Nk     3#idin  4#-1      5#+1    6aIDin        7IDall         8FS+1          9FS-1          10FS           11ideal       12eidin  13eidall" > $af
fi
awk '{print $0}' fort.35 >> $af
### @@@@-U
#obtains number of k-points
if [ ! -f fort.246 ]
   then
       sh linea.sh
       sh warning.sh
       printf "\tIncrease the value of -k so there k-points\n"
       sh linea.sh
       exit 1
fi
sh linea.sh
#converts cartesian to direct k-points of the Fermi Surface
if [ "1" == "1" ]
then
    sh linea.sh
    printf "\tOutput: Original $N-divisions  => the following Cartesian IBZ-kpoints\n"
    if (( $(echo "$kcBand > .0" | bc -l) ))
    then #band 36 or 37
	printf "\tfor Band->$band \n"
	##
	#gets info
	Ncg=`awk '{print $1}' .info`
	Nr=`awk '{print $2}' .info`
	#reading for the 2nd line
	kni=`awk 'FNR==2{print $2}' .coarse-grained_$Nr`
	#flash.sh "AquiBoyNas"
	#printf "\tNcg=$Ncg Nr=$Nr kni=$kni\n"
	cual=$sym/$case.fs-ibz-in-cg-$Nr-r-$Nr-k-$kni
	echo $Ncr $Nr $kni >> .info-cg-r-kf
	#ene=`wc fort.246 | awk '{print $1}'` #number of FS k-points
	#cual=$case.kflist_$ene
	#cp fort.246 $cual
	##
	printf "\t@$script: Converting Cartesian all-kpoints of $cual to direct lattice kpoints\n"
	#
	awk '{print $1,$2,$3}' $cual > .auxfs
	#Converts cartesian to direct k-points
	#kcart2dlatt.sh -f k-file -o (1,2)[Cartesian->DL,DL->Cartesian]
	sh linea.sh
	#printf "\tkcart2dlatt.sh -f .auxfs -o 1\n"
	kcart2dlatt.sh -f .auxfs -o 1
	Nkfs=`wc .auxfs-direct | awk '{print $1}'`
#	echo $Nkfs > .knfs_$nk
	out21d=$case.klist_$Nkfs
	mv .auxfs-direct $out21d
    else #band 34 or 35
	printf "\tfor Band->$band\n"
	##
	#gets info
	Ncg=`awk '{print $1}' .info`
	Nr=`awk '{print $2}' .info`
	kni=`awk '{print $2}' fort.123`
	cual=$sym/$case.fs-ibz-in-cg-$Nr-r-$Nr-k-$kni
	echo $Ncg $Nr $kni >> .info-cg-r-kf
	##
	printf "\t@$script: Converting Cartesian k-points of $cual to direct lattice kpoints\n"
	#
	awk '{print $1,$2,$3}' $cual > .auxfs
	#Converts cartesian to direct k-points
	#kcart2dlatt.sh -f k-file -o (1,2)[Cartesian->DL,DL->Cartesian]
	sh linea.sh
	kcart2dlatt.sh -f .auxfs -o 1
	Nkfs=`wc .auxfs-direct | awk '{print $1}'`
#	echo $Nkfs > .knfs_$nk
	out21d=$case.klist_$Nkfs
	mv .auxfs-direct $out21d
    fi
    #
    flash.sh "resulting in:"
    sh linea.sh
    printf "\t${blue}$out21d${NC}  -> ${red}Direct IBZ $Nkfs k-points to be used by run_tiniba.sh for the Fermi-Surface${NC}\n"
    knmin=`awk '{print $1}' fort.248`
    knmax=`awk '{print $2}' fort.248`
    knmin2=`awk '{print $1}' fort.247`
    knmax2=`awk '{print $2}' fort.247`
    if [[ "$band" == "34" || "$band" == 35 ]]
       then
	   printf "\t${blue}$finm${NC} -> ${red}$nm Cartesian IBZ Coarse-Grained k-points${NC}\n"
	   printf "\t${CYAN}$cual${NC} -> ${red}$Nkfs Cartesian IBZ next Refined k-points${NC}\n"
	   printf "\t${CYAN}.coarse-grained_$Nr${NC} -> ${red}with info for next Refined k-points${NC}\n"
    fi
    if [[ "$band" == "36" || "$band" == 37 ]]
       then
	   printf "\t${blue}$finm${NC} -> ${red}$Ncg Cartesian IBZ Coarse-Grained k-points${NC}\n"
	   printf "\t                                                   ${red}with $nm grad=-1${NC}\n"
	   printf "\t${blue}$finp${NC} -> ${red}$Ncg Cartesian IBZ Coarse-Grained k-points${NC}\n"
	   printf "\t                                                  ${red}with $np grad=+1${NC}\n"
	   printf "\t${CYAN}$cual${NC} -> ${red}$Nkfs Cartesian IBZ next Refined k-points${NC}\n"
	   printf "\t${CYAN}.coarse-grained_$Nr${NC} -> ${red}with info for next Refined k-points${NC}\n"
	   #printf "\t                   kFS_min-1=$knmin & kFS_max-1=$knmax\n"
	   #printf "\t                   kFS_min-2=$knmin2 & kFS_max-2=$knmax2\n"
    fi
    cuale=eigen_$Nkfs\_$ecut-$spin
    sh linea.sh
    printf "\t${CYAN}Calculating $cuale & corresponding momentum matrix elements${NC}\n"
    sh linea.sh
    ######################################################
    #################Exit-Down
    #################################
    if [ "1" == "2" ]               #
    then                            #
	sh linea.sh                 #
	flash.sh "@$script:Exit#7"  #
	sh linea.sh                 #
	exit 1                      #
    fi                              #
    #################################
    #################Exit-Up
    ################################# run_tiniba.sh-Down #
    if [ ! -f $cuale ]
    then
	#### Commands
	nodo=`hostname`
	if [[ "$nodo" == "fat1" || "$nodo" == "fat2" || "$nodo" == "fat3" ]]
	then
	    cores=64
	else
	    cores=12
	fi
	if [[ "$nodo" == "viginti01" || "$nodo" == "viginti02" ]]
	then
	    cores=40
	fi
	#Number of plane waves
	pw=`grep mpw $case'_check'/$case.out | awk '{print $12}'`
	sh linea.sh
	printf "\t@$script: using $cores cores for SCF\n"
	sh linea.sh
	#here we  need -p and -v, since these are the refined FS-k-points
	exectiniba="run_tiniba.sh -r run -k $Nkfs  -N 0 -x 2 -C $cores -P $pw -w -e -p -v"
	printf "\t${blue}####### executing for energies and momentum matrix elements (including Vnl) for $Nkfs FS-kpoints ######${NC}\n"
	printf "\t${CYAN}$exectiniba\n"
	printf "\t${blue}#####################################################################################${NC}\n"
	# executable
	$exectiniba
	#####-Down
	# this will be deprecated once we find how to compile DP->Cabellos work!
	sh linea.sh
	exectiniba="run_tiniba.sh -r run -k $Nkfs  -N 0 -x 2 -C $cores -P $pw -Z"
	printf "\t${blue}####### executing for diagonal momentum matrix elements (including Vnl) for $Nkfs FS-kpoints ######${NC}\n"
	printf "\t${CYAN}$exectiniba\n"
	printf "\t${blue}#####################################################################################${NC}\n"
	# executable
	#################Exit-Down
	if [ "1" == "2" ]
	then
	    sh linea.sh
	    flash.sh "@$script:Exit#8"
	    sh linea.sh
	    exit 1
	fi
	#################Exit-Up
	$exectiniba
	#####-UP
	sh linea.sh
	printf "\t${cyan}run_tiniba.sh done!${NC} back @$script\n"
	sh linea.sh
    else
	flash.sh "case $cuale done!"
	sh linea.sh
    fi
    ################################# run_tiniba.sh-Up ##
    #####################################################
    sh linea.sh
    printf "\tPlot:\n"
    printf "\tFull BZ\n"
    printf "\tgnuplot> load 'symmetries/bz-bcc.g'\n"
    printf "\t${RED}or${NC}\n"
    printf "\tIBZ\n"
    printf "\tgnuplot> load 'symmetries/ibz-bcc.g'\n\n"
    printf "\t${blue}Fermi Surface${NC}\n"
    # all-k-points for FS
    if (( $(echo "$kcBand > .0" | bc -l) ))
    then
	printf "\tgnuplot> sp '$cual' u (\$9==-1?\$1:1/0):2:3 w p pt 7 ps .9 lc rgb \"blue\"  t 'Refined:-1'\n"
	printf "\tgnuplot> rep '$cual' u (\$9==1?\$1:1/0):2:3 w p pt 7 ps .9 lc rgb \"red\"   t 'Refined:+1'\n"
	printf "\tgnuplot> rep '$finm' u 1:2:3 w p pt 7 ps .9 lc rgb \"green\"  t 'CG:-1'\n"
	printf "\tgnuplot> rep '$finp' u 1:2:3 w p pt 7 ps .9 lc rgb \"cyan\"  t 'CG:+1'\n"
	printf "\tDuck@set view 90,225,2.5\n"
	printf "\tOriginal@set view 82,119,2.5\n"
	# Octant Plots
	sh linea.sh
	printf "\tOctant Plots:\n"
	goutf=$outdir/octant-fs-$N-$Ncg.g
	printf "\tgnuplot> load 'symmetries/bz-bcc.g'\n"
	printf "\t${RED}or${NC}\n"
	printf "\tgnuplot> load 'symmetries/ibz-bcc.g'\n\n"
	printf "\tgnuplot> load '$goutf'\n"
	echo "set ticslevel 0" > $goutf
	sp=.7
	cl=1
	echo "sp '$cual' u (\$9==-1?\$1:1/0):3:2 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "   '$cual' u (\$9==-1?\$2:1/0):3:1 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "   '$cual' u (\$9==-1?\$2:1/0):1:3 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "   '$cual' u (\$9==-1?\$3:1/0):1:2 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "   '$cual' u (\$9==-1?\$3:1/0):2:1 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "   '$cual' u (\$9==-1?\$1:1/0):2:3 w p pt 7 ps $sp lc $cl t 'Refined:-1',\\" >> $goutf
	cl=2
	echo "   '$cual' u (\$9==1?\$1:1/0):3:2 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "   '$cual' u (\$9==1?\$2:1/0):3:1 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "   '$cual' u (\$9==1?\$2:1/0):1:3 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "   '$cual' u (\$9==1?\$3:1/0):1:2 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "   '$cual' u (\$9==1?\$3:1/0):2:1 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "   '$cual' u (\$9==1?\$1:1/0):2:3 w p pt 7 ps $sp lc $cl t 'Refined:+1',\\" >> $goutf
	cl=3
	echo "'$finm' u 1:3:2 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "'$finm' u 2:3:1 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "'$finm' u 2:1:3 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "'$finm' u 3:1:2 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "'$finm' u 3:2:1 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "'$finm' u 1:2:3 w p pt 7 ps $sp lc $cl t 'CG:-1',\\" >> $goutf
	cl=4
	echo "'$finp' u 1:3:2 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "'$finp' u 2:3:1 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "'$finp' u 2:1:3 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "'$finp' u 3:1:2 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "'$finp' u 3:2:1 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "'$finp' u 1:2:3 w p pt 7 ps $sp lc $cl t 'CG:+1'" >> $goutf
    else
	printf "\tgnuplot> sp '$out4' u (\$6==0?\$1:1/0):2:3 w p pt 7 ps .6  t 'BF=0'\n"
	printf "\tgnuplot> rep '$out4' u (\$6==1?\$1:1/0):2:3 w p pt 7 ps .6  t 'BF=1'\n"
	printf "\tgnuplot> sp '$finm' u 1:2:3 w p pt 7 ps .9 lc rgb \"blue\"  t 'CG-FS'\n"
	printf "\tgnuplot> rep '$cual' u 1:2:3 w p pt 7 ps .9 lc rgb \"red\"   t 'Refined-FS'\n"
	# Octant Plots
	sh linea.sh
	printf "\tOctant Plots:\n"
	goutf=$outdir/octant-fs-$N-$Ncg.g
	printf "\tgnuplot> load 'symmetries/bz-bcc.g'\n"
	printf "\t${RED}or${NC}\n"
	printf "\tgnuplot> load 'symmetries/ibz-bcc.g'\n\n"
	printf "\tgnuplot> load '$goutf'\n"
	echo "set ticslevel 0" > $goutf
	cl=1
	sp=.7
	echo "sp '$cual' u 1:3:2 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "   '$cual' u 2:3:1 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "   '$cual' u 2:1:3 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "   '$cual' u 3:1:2 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "   '$cual' u 3:2:1 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "   '$cual' u 1:2:3 w p pt 7 ps $sp lc $cl t 'Refined-FS',\\" >> $goutf
	cl=2
	echo "'$finm' u 1:3:2 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "'$finm' u 2:3:1 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "'$finm' u 2:1:3 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "'$finm' u 3:1:2 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "'$finm' u 3:2:1 w p pt 7 ps $sp lc $cl t '',\\" >> $goutf
	echo "'$finm' u 1:2:3 w p pt 7 ps $sp lc $cl t 'CG-$Ncg'" >> $goutf
    fi
    sh linea.sh
fi
sh linea.sh
flash.sh "Data output@:"
printf "\t$af\n"
sh linea.sh
sh end.sh
sh linea.sh
printf "\t@$script: $como_corri\n"
sh linea.sh
#rm -f fort.* sort* dog .info .knfs_ sera.dat
rm -f fort.* sort* dog .info sera.dat
exit
#########################################
minimorum-401

cp ../../eigen_1906_14-nospin .
cp ../../.coarse-grained_351 .
cp ../../setUpAbinit_y2c3.in .
cp ../../symmetries/y2c3.fs-ibz-in-cg-351-r-351-k-345 symmetries/.
cp ../../.spininfo .

401/copy-401.sh
451/copy-451.sh




##################################################################################################
band 34
symmetries/y2c3.fs-ibz-in-179  -> Direct IBZ 179 k-points to be used by run_tiniba for the Fermi-Surface
	                          kFS_min=0.040993 & kFS_max=0.046323

band35
symmetries/y2c3.fs-ibz-in-172  -> Direct IBZ 172 k-points to be used by run_tiniba for the Fermi-Surface
	                          kFS_min=0.046245 & kFS_max=0.061085

band 36
symmetries/y2c3.fs-ibz-in-108  -> Direct IBZ 108 k-points to be used by run_tiniba for the Fermi-Surface
	                          kFS_min-1=0.158662 & kFS_max-1=0.200421
	                          kFS_min-2=0.064607 & kFS_max-2=0.085186

band 37
symmetries/y2c3.fs-ibz-in-364  -> Direct IBZ 364 k-points to be used by run_tiniba for the Fermi-Surface
	                          kFS_min-1=0.131189 & kFS_max-1=0.236775
	                          kFS_min-2=0.067998 & kFS_max-2=0.104791
