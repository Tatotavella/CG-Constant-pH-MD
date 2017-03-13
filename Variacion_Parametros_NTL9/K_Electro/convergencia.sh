#!/bin/bash

# Analisis de la convergencia de la fraccion desprotonada
LANG=en_US

freqMC=5


totalResCharged=15
columns=$((totalResCharged+1))
in=ntl9
MCfile=MC.input
MCres=MC.result
MCdeprot=MC.deprot.ejemplo.4.15
out=convergencia
pHini=1
pHfin=14
pasopH=1


rm $out.txt

#res=`awk 'NR==1{print}' $MCdeprot`
#echo -e $res >> "$out".txt 

while read -r var
do
write=""
rm uso.txt
echo -e $var >> uso.txt 
for j in $(seq 1 1 $columns); do
  res=`awk 'NR==1{print}' uso.txt | awk -v var=$j '{print $var/($1+1)}'`
  #echo $res
  write+="|$res"
done
echo -e $write >> "$out".txt 
done < "$MCdeprot"

column -s"|" -t "$out".txt >> salida.txt
mv salida.txt "$out".txt





#write="Step\Resid"
#for j in $(seq 1 1 $totalResCharged); do
#res=`awk 'NR==1{print}' $MCdeprot`
# | awk -v var=$j '{print $var/($1+1)}'`
#echo $res
#write+="|$res"
#done

#echo -e $write >> "$out".txt 

#lineEnds=""

