
model Basic -ndm 3 -ndf 6

node 1  0.0               0.0  0.0 
node 2  9.801714032956060 0.0  0.4815273327803071 
node 3 19.509032201612835 0.0  1.9214719596769570
node 4 29.028467725446234 0.0  4.3059664267791180
node 5 38.268343236508976 0.0  7.6120467488713260
node 6 47.139673682599760 0.0 11.8078735651644950
node 7 55.557023301960214 0.0 16.8530387697454780
node 8 63.439328416364546 0.0 22.6989546637263000
node 9 70.710678118654740 0.0 29.2893218813452400
fix 1 1 1 1 1 1 1


geomTransf Corotational 1 0.0 0.0 1.0 

set E   1000.0
set A  10000.0
set Ay 10000.0
set Az 10000.0
set I    833.3333333333334
set G 500.0
set J 1666.6666666666667

element elasticBeamColumn 1 1 2 $A $E $G $J $I $I 1
element elasticBeamColumn 2 2 3 $A $E $G $J $I $I 1
element elasticBeamColumn 3 3 4 $A $E $G $J $I $I 1
element elasticBeamColumn 4 4 5 $A $E $G $J $I $I 1
element elasticBeamColumn 5 5 6 $A $E $G $J $I $I 1
element elasticBeamColumn 6 6 7 $A $E $G $J $I $I 1
element elasticBeamColumn 7 7 8 $A $E $G $J $I $I 1
element elasticBeamColumn 8 8 9 $A $E $G $J $I $I 1


pattern Plain 1 Linear {
  load 9  0.0 0.0 0.0 0 0 0 ;
}

system BandGeneral 
numberer RCM 
constraints Plain 
algorithm Newton 
test NormUnbalance 1e-6 20 

integrator LoadControl 0.0
analysis Static 
analyze 1 

puts "  Iterations: [numIter]"

for {set i 1} {$i < 10} {incr i} {
  puts "  Displacement: [nodeDisp $i 1] [nodeDisp $i 2] [nodeDisp $i 3]" 
}

