# OpenSees Examples

This directory contains 5 examples that demonstrate improvements of the proposed transformations over the existing `CrdTransf` classes in OpenSees.

All tests can be run in a Posix shell by sourcing the script `test.sh`.

<dl>
<dt><b>Test 01</b></dt>
<dd>
This example demonstrates the convergence characteristics of the new corotational transformations.
The standard problem of a curved 45-degree cantilever is implemented.
The following variants are investigated:
 <dl>
 <dt><tt>Test01-Geom01</tt></dt>
 <dd>This is the formulation that is currently available in OpenSees as the <tt>Corotational</tt> transformation.
     It can be executed with both the <tt>OpenSees</tt> executable, and the <tt>xara</tt> executable.
 </dd>
 <dt><tt>Test01-Geom02</tt></dt>
 <dd>
 </dd>
 </dl>
</dd>
<dt><b>Test 06</b></dt>
<dd>This example demonstrates the use of the Corotational02 transformation to represent a shear-deformable cantilever.
The setup is that of Section 4.2.2 from Perez and Filippou (2024).
 <dl>
 <dt><tt>Test01-Prism-Geom01</tt></dt>
 <dd>This is the formulation that is currently available in OpenSees as the <tt>Corotational</tt> transformation.
     It can be executed with both the <tt>OpenSees</tt> executable, and the <tt>xara</tt> executable.
 </dd>
 <dt><tt>Test01-Prism-Geom02</tt></dt>
 <dd>
 </dd>
 </dl>
</dd>
</dl>
