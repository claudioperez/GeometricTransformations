
model Basic -ndm 3 -ndf 6
set E 10e3
set G 10e3
set A 1
set I 1e-2
set J 1e-2

section ElasticFrame 1 $E $A $I $I $G $J 1.0 1.0

geomTransf Corotational02 1 0 0 1 

node  0  0.0 0 0 
node  1  1.0 0 0 
node  2  2.0 0 0 
node  3  3.0 0 0 
node  4  4.0 0 0 
node  5  5.0 0 0 
node  6  6.0 0 0 
node  7  7.0 0 0 
node  8  8.0 0 0 
node  9  9.0 0 0 
node 10 10.0 0 0 

element PrismFrame  1 0  1  -section 1 -transform 1
element PrismFrame  2 1  2  -section 1 -transform 1
element PrismFrame  3 2  3  -section 1 -transform 1
element PrismFrame  4 3  4  -section 1 -transform 1
element PrismFrame  5 4  5  -section 1 -transform 1
element PrismFrame  6 5  6  -section 1 -transform 1
element PrismFrame  7 6  7  -section 1 -transform 1
element PrismFrame  8 7  8  -section 1 -transform 1
element PrismFrame  9 8  9  -section 1 -transform 1
element PrismFrame 10 9 10  -section 1 -transform 1
fix  0 1 1 1 1 1 1 
fix 10 0 0 0 0 0 0 


pattern Plain 1 Linear {
  load 10 0 0 25 0 0 314.1592653589793
}

set nstep 400
system Umfpack 
test NormUnbalance 1e-10 55 0 
numberer RCM 
constraints Plain 
algorithm Newton 
integrator LoadControl [expr 1.0/$nstep] 
analysis Static 

analyze $nstep
puts "  Iterations: [numIter]"
puts "  Time: [getTime]"
puts "  Displacement: [nodeDisp 10 1] [nodeDisp 10 2]"

