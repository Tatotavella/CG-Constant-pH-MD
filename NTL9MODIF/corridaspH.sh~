#!/bin/bash
LANG=en_US

freqMC=5

pKa1=4.5
pKa2=5.0

totalResCharged=15
columns=$((totalResCharged+1))
in=testeo
MCfile=MC.input
MCres=MC.result
MCdeprot=MC.deprot
out=result


rm $out.txt

write="pH"
for j in $(seq 1 1 $totalResCharged); do
res=`awk 'NR>3{print $1}' $MCfile | awk -v var=$j 'FNR==var{print}'`
echo $res
write+="|$res"
done

echo -e $write >> "$out".txt 



for i in $(seq 2 2 12); do
write=""
write+="$i"
sed "s/freqMC/$freqMC/g" "$MCfile".bk | sed "s/pH/$i/g"  > $MCfile
#sed "s/PH/$i/g" "$in".in.bk > "$in".in
./lmp_serial < "$in".in
endofline=`awk 'END{print}' $MCdeprot`
echo $endofline
for j in $(seq 2 1 $columns); do
deprot=`awk 'END{print}' $MCdeprot | awk -v var=$j '{print $var/($1+1)}'`
#echo $deprot
write+="|$deprot"

echo -e $write >> "$out".txt 
column -s"|" -t "$out".txt >> salida.txt
done
mv salida.txt "$out".txt






#total1=`awk -v var=$r1 'NR>1{if ($1==var) print $2,$3}' $MCdeprot | wc | awk '{print $1}'`


#deprot1=`awk -v var=$r1 -v vare=$r1q 'NR>1{if ($1==var && $2==vare) print $2,$3}' $MCres | wc | awk '{print $1}'`
#total2=`awk -v var=$r2 'NR>1{if ($1==var) print $2,$3}' $MCres | wc | awk '{print $1}'`
#deprot2=`awk -v var=$r2 -v vare=$r2q 'NR>1{if ($1==var && $2==vare) print $2,$3}' $MCres | wc | awk '{print $1}'`
#echo fd=$((deprot/total))
#echo $i $deprot1 >> "$r1".txt
#echo $i $deprot2 >> "$r2".txt
done
