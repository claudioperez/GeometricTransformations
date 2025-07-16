
xara=/Users/claudio/packages/OpenSeesRT/build/temp.macosx-11.0-arm64-cpython-312_local/SRC/executable/xara
OpenSees=/Users/claudio/opensees/OpenSees/build-py3132/OpenSees #/Users/claudio/opensees/OpenSees/build/OpenSees
hyperfine --warmup 20 --runs 50 \
 "              $OpenSees Test05-Prism-Geom01.tcl -noHeader" \
 "              $xara     Test05-Prism-Geom01.tcl" \
 "              $xara     Test05-Prism-Geom02.tcl" \
 "Crisfield=1   $xara     Test05-Prism-Geom02.tcl" \
 "Crisfield02=1 $xara     Test05-Prism-Geom02.tcl"

