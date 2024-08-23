# `GeometricTransformations`


A framework for geometrically exact finite element analysis

[![DOI:10.1002/nme.7506](https://zenodo.org/badge/DOI/10.1002/nme.7506.svg)](https://doi.org/10.1002/nme.7506)

--------------------------------------------------

> **This is a partial repository; an upload of the full project is in progress**

This project implements a framework for simulating differential equations on
nonlinear manifolds. 
<!--
Theoretical developments give rise to a modular computational
framework for composing coordinate transformations and manifold
parameterizations. 
-->
This repository contains the source code for the investigation presented by Perez and Filippou (2024).
The MATLAB implementation requires the [FEDEASLab](https://fedeas.net) toolbox for nonlinear structural analysis
which may be obtained from the author upon reasonable request.
A completely open-source C++ implementation of the framework is under development 
for OpenSees and expected to be completed by 2025.

<dl>
<dt><a href="./Elements">Elements/</a></dt>
<dd>This directory contains the finite elements:
  <ul> 
    <li><a href="./Elements/DisplShear3dFrm_wCS.m"><code>DisplShear3dFrm_wCS</code></a> 
       Displacement-interpolated shear 3d frame element with Cosserat strains; 
       this element implements the three geometrically exact beam formulations 
       presented in Appendix B.
    </li>
    <li><a href="./Elements/GeomWrap3dFrm.m"><code>GeomWrap3dFrm</code></a> Wrapper for 3d frame elements. This element
       implements the <em>element wrapper</em> from Section 6 using Algorithm 1
       from Section 4.
    </li>
    <li><a href="./Elements/GeomTran3dFrm.m"><code>GeomTran3dFrm</code></a> Wrapper for 3d frame elements
       that are formulated in a <em>basic</em> coordinate system. This element
       implements the <em>element wrapper</em> from Section 6 using Algorithm 1
       from Section 4.
    </li>
  </ul>
</dd>
<dt><a href="./Geometry">Geometry/</a></dt>
<dd>
  <ul>
    <li><a href="/claudioperez/FiniteRotationLab/blob/master/Geometry/GeomTran3dFrm_Pull.m"><code>GeomTran3dFrm_Pull</code></a>
    This function implements the $g$ operation of the geometric transformations by calling the following functions:
    <ul>
      <li><a href="/claudioperez/FiniteRotationLab/blob/master/Geometry/Transform_RotPull.m"><code>Transform_RotPull</code></a></li>
      <li><a href="/claudioperez/FiniteRotationLab/blob/master/Geometry/Transform_IsoPull.m"><code>Transform_IsoPull</code></a></li>
      <li><a href="/claudioperez/FiniteRotationLab/blob/master/Geometry/Tran3dFrm_IsoPull.m"><code>Tran3dFrm_IsoPull</code></a></li>
    </ul>
    </li>
    <li><a href="/claudioperez/FiniteRotationLab/blob/master/Geometry/GeomTran3dFrm_Push.m"><code>GeomTran3dFrm_Push</code></a>
    This function implements the $\mathbf{a}_g$ and $\mathbf{k}_g$ operations of the geometric transformations.
    <ul>
      <li><a href="/claudioperez/FiniteRotationLab/blob/master/Geometry/Tran3dFrm_IsoPush.m"><code>Tran3dFrm_IsoPush</code></a></li>
    </ul>
    </li>
  </ul>
</dd>
<dt><a href="./Rotations">Rotations/</a></dt>
<dd>This directory contains the rotation functions described in Appendix A.
</dd>

<dt><a href="./Examples">Examples/</a></dt>
<dd>
  This directory contains scripts to reproduce the examples of the paper.
  <ol>
  <li><a href="./Examples/E10_Invariance.m">Invariance</a></li>
  <li><ol>
      <li><a href="./Examples/E21_PlaneMoment.m">Plane Cantilever - Point moment</a></li>
      <li><a href="./Examples/E22_PlaneTransverse.m">Plane Cantilever - Transverse Force</a></li>
  </ol></li>
  <li><a href="./Examples/E30_HelicalForms.m">Space Cantilever - Point moment and force</a></li>
  <li><a href="./Examples/E40_BatheCantilever.m">Curved Cantilever - Point force (Bathe's problem)</a></li>
  <li><a href="./Examples/E50_Hockling.m">Hockling of a flexible rod.</a></li>
  </ol>
</dd>
</dl>

For an example using "basic" elements, see <a href="./Examples/E00_Cantilever.m"><code>E00_Cantilever</code></a>.

<table align="center" style="border: 0;">
 <tr style="background-color:rgba(0, 0, 0, 0);">
  <td style="background-color:rgba(0, 0, 0, 0);" colspan="3">
    <a>
    <img src="./Figures/Figure_4a.png" 
         width="600" alt="OpenSeesRT Logo">
    </a>
  </td>
 </tr>
</table>


- Perez CM, Filippou FC. On nonlinear geometric transformations of finite elements. Int
  J Numer Methods Eng. 2024;e7506. doi: 10.1002/nme.7506

