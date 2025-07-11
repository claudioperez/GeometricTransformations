#alias OpenSees="OpenSees -noHeader"


for i in 1 2 3 4; do
  cat <<END

----------------
Test $i
----------------
END
  printf "\na) Corotational (OpenSees)\n"
  OpenSees 0$i/Test0$i-Force-Geom01.tcl -noHeader 
  printf "\nb) Corotational   (xara)\n"
  xara     0$i/Test0$i-Force-Geom01.tcl
  printf "\nc) Corotational02 (xara)\n"
  xara     0$i/Test0$i-Force-Geom02.tcl
  printf "\nc) Corotational03 (xara)\n"
  Crisfield=1 xara 0$i/Test0$i-Force-Geom02.tcl
done


cat <<END

----------------
Test 5
----------------
END
time (repeat 10 { OpenSees Test05-Force-Geom01.tcl -noHeader 2>/dev/null; } )
time (repeat 10 { xara     Test05-Force-Geom01.tcl; } )
time (repeat 10 { xara     Test05-Force-Geom02.tcl; } )
time (repeat 10 { Crisfield=1 xara     Test05-Force-Geom02.tcl; } )

