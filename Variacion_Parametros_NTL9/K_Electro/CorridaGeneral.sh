#!/bin/bash

LANG=en_US



lammps=/home/er17/lammps-9Oct12/src


#Archivos a modificar
in=ntl9
backbone=fix_backbone_coeff



dir=`pwd`

#Modificaciones de ntl9.in
TSTEP=2
THERMO=10
RUNTIME=50000
DUMP=50
TINI=300.0
TFIN=300.0


#Modificaciones de fix_backbone_coeff.data

#Constante de Debye
CONTACTOINI=4.15 #Este es el que va a cambiar
CONTACTOSTEP=1
CONTACTOFIN=5.15

#Apantallamiento Debye
PANTALLA=10.0


#Para cada valor de la constante electrostatica
#for C in $(seq $CONTACTOINI $CONTACTOSTEP $CONTACTOFIN); do
for C in 0 1;do
	rm -r $C
	mkdir $C
	cd $C
	cp $dir/files/* .
	SEMILLA=`shuf -i1000000-9999999 -n1` #Numero al azar de 7 cifras
	#Busco las palabras para cambiar
	sed "s/XXXXXXX/$SEMILLA/g" "$in"bak.in | sed "s/TSTEP/$TSTEP/g" | sed "s/DUMP/$DUMP/g" | sed "s/THERMO/$THERMO/g" | sed "s/RUNTIME/$RUNTIME/g" | sed "s/TFIN/$TFIN/g" | sed "s/TINI/$TINI/g" > "$in".in
	sed "s/PANTALLA/$PANTALLA/g" "$backbone"bak.data | sed "s/CONTACTO/$C/g" > "$backbone".data
	cd ..
	pwd
done
pwd

#for C in $(seq $CONTACTOINI $CONTACTOSTEP $CONTACTOFIN); do echo cd $C; echo cp $dir/files/run.sh $lammps/lmp_serial .;echo "sbatch run.sh >> $dir/submitted"; echo "pwd >> $dir/largadas";echo cd ..; done; echo cd ..; done; echo cd ..; done >> biou.sh


#chmod +x biou.sh

#echo ./biou.sh





