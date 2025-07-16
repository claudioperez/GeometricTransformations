# Pure flexure of a cantilever
set element forceBeamColumn

model Basic -ndm 3 -ndf 6

node 1 0.0 0.0 0.0 
node 2 0.2 0.0 0.0 
node 3 0.4 0.0 0.0 
node 4 0.6 0.0 0.0 
node 5 0.8 0.0 0.0 
node 6 1.0 0.0 0.0 

fix 1 1 1 1 1 1 1 
fix 6 0 0 1 1 1 0 

set E  1.0
set G  1.0 
set A  2.0
set J  2.0
set Iy 2.0
set Iz 2.0
set Ay 2.0
set Az 2.0

section Elastic 1 $E $A $Iy  $Iz $G $J 1.0 1.0

geomTransf Corotational 1 0 0 1 

element $element 1  1 2  5  1 1
element $element 2  2 3  5  1 1
element $element 3  3 4  5  1 1
element $element 4  4 5  5  1 1
element $element 5  5 6  5  1 1


test EnergyIncr 1e-12 10 0
pattern Plain 1 Linear {
  load 6 0 0 0 0 0 [expr -4.0*acos(-1)];
}

integrator LoadControl 0.2
constraints Plain
numberer RCM
system BandGen
algorithm Newton
analysis  Static 
analyze 5

puts "Iterations: [numIter]"
puts "Displacement: [nodeDisp 6 1]"
