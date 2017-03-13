#!/bin/bash
LANG=en_US

freqMC=5


totalResCharged=15
columns=$((totalResCharged+1))
in=ntl9
MCfile=MC.input
MCres=MC.result
MCdeprot=MC.deprot
out=result
pHini=1
pHfin=14
pasopH=1


rm $out.txt

write="pH"
for j in $(seq 1 1 $totalResCharged); do
res=`awk 'NR>3{print $1}' $MCfile | awk -v var=$j 'FNR==var{print}'`
echo $res
write+="|$res"
done

echo -e $write >> "$out".txt 

lineEnds=""

for i in $(seq $pHini $pasopH $pHfin); do
write=""
write+="$i"
sed "s/freqMC/$freqMC/g" "$MCfile".bk | sed "s/pH/$i/g"  > $MCfile
#sed "s/PH/$i/g" "$in".in.bk > "$in".in
./lmp_serial < "$in".in


endofline=`awk 'END{print}' $MCdeprot`
lineEnds+="\n$endofline"


for j in $(seq 2 1 $columns); do
deprot=`awk 'END{print}' $MCdeprot | awk -v var=$j '{print $var/($1+1)}'`
write+="|$deprot"
done

echo $write >> "$out".txt 

done
echo -e $lineEnds
column -s"|" -t "$out".txt >> salida.txt
mv salida.txt "$out".txt
