#!/bin/bash
LANG=en_US

freqMC=5

pKa1=4.5
pKa2=5.0
in=testeo
MCfile=MC.data
r1=22
r2=37
r1q=-1
r2q=-1
MCres=MC.result
out=res.txt
rm $out $r2.txt $r1.txt

for i in $(seq 2 2 12); do

sed "s/pKa1/$pKa1/g" "$MCfile".bk | sed "s/pKa2/$pKa2/g" | sed "s/pH/$i/g" | sed "s/freqMC/$freqMC/g" > $MCfile
sed "s/PH/$i/g" "$in".in.bk > "$in".in
./lmp_serial < "$in".in
total1=`awk -v var=$r1 'NR>1{if ($1==var) print $2,$3}' $MCres | wc | awk '{print $1}'`
deprot1=`awk -v var=$r1 -v vare=$r1q 'NR>1{if ($1==var && $2==vare) print $2,$3}' $MCres | wc | awk '{print $1}'`
total2=`awk -v var=$r2 'NR>1{if ($1==var) print $2,$3}' $MCres | wc | awk '{print $1}'`
deprot2=`awk -v var=$r2 -v vare=$r2q 'NR>1{if ($1==var && $2==vare) print $2,$3}' $MCres | wc | awk '{print $1}'`
#echo fd=$((deprot/total))
echo $i $deprot1 >> "$r1".txt
echo $i $deprot2 >> "$r2".txt
done
