/** @file tutorial.hh 
 * @todo separate into multiple files named *.dox
 */
/**
\mainpage

The Sparse Least squares on Manifolds (SLoM) package aims to make optimizing 
least squares problems arising, e.g., in SLAM or calibration problems straight
forward and reasonably fast. It furthermore exploits the @ref MANIFOLD "manifold" 
structure of variables such as orientations (@f$\SO n@f$) or unit vectors 
(@f$S^2@f$), thus elegantly avoiding problems with singularities or over-parameterization.


\section REQ Requirements
The Manifold ToolKit requires the Boost/Preprocessor library and the Eigen Matrix library 
(version 2.0.9 or later &mdash; Eigen 3 is not supported at the moment).

SLoM furthermore requires the CXSparse library (often available in a package called 'SuiteSparse'),
Boost/Intrusive and Boost/Bind.

SLoM and MTK have been developed and tested with g++ (3.4.4 and later), if you have problems or 
need help porting to other compilers, please contact the author.


\section Building
SLoM and MTK are header-only libraries, so no building is required.
There is a CMakeFile to build the examples and unit-tests (see @ref BUILD "here" for details).
Compiling own programs requires MTK, SLoM and the @ref REQ "required libraries" in the include path,
linking needs -lcxsparse (and corresponding library in link path).


\section FS First steps

First steps to understand SLoM and MTK is to read the @ref SLOMTutorial, for more details read the
@ref MTKTutorial.


\section References

See \ref bibpage for a list of citations. When using SLoM, please cite one of the following:
Further (more theoretical) details are provided in:
\verbatim
@MastersThesis{Hertzberg2008,
  author = {C. Hertzberg},
  title  = {A Framework for Sparse, Non-Linear Least Squares Problems on Manifolds},
  school = {Universit\"at Bremen},
  year   = {2008},
}

@Article{Hertzberg2010,
  author  = {C. Hertzberg and R. Wagner and U. Frese and L. Schr\"oder},
  title   = {Integrating Generic Sensor Fusion Algorithms with 
             Sound State Representations through Encapsulation of Manifolds},
  journal = {Information Fusion},
  year    = {submitted March 2010},
  note    = {\textbf{to appear}},
}
\endverbatim


\section License

 *  Copyright (c) 2008&ndash;2011, Universit&auml;t Bremen
 *  All rights reserved.
 *
 *  Authors: 
 *   - Christoph Hertzberg <chtz@informatik.uni-bremen.de> and
 *   - Ren&eacute; Wagner <rene.wagner@dfki.de>
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   - Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   - Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   - Neither the name of the Universitaet Bremen nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.

*/

/**
\page bibpage Bibliography

@biblio
 @bibentry{Hertzberg2008, C. Hertzberg, A Framework for Sparse\, Non-Linear Least Squares Problems on Manifolds, Diplomarbeit, Univ. Bremen 2008}
 @bibentry{Hertzberg2011, C. Hertzberg\, R. Wagner\, U. Frese\, L. Schr&ouml;der, Integrating Generic Sensor Fusion Algorithms with Sound State Representations through Encapsulation of Manifolds, Information Fusion, to appear}
 @bibentry{Lee2003, J.M. Lee, Introduction to Smooth Manifolds, volume 218 of Graduate Texts in Mathematics, Springer Verlag 2003}
@endbiblio



*/


/**
\page SLOMTutorial Tutorial to SLoM

This page shows a minimal example how to initialize and optimize a pure 
pose relation problem.

\include example.cpp
*/


/**
\page MTKTutorial Tutorial to MTK

The Manifold ToolKit (MTK) provides convenient ways to handle @ref MANIFOLD "manifolds"
as they occur in many estimation or sensor fusion problems.
At first MTK provides some manifold primitives such as plain vectors @f$\R^n@f$ 
(see @c MTK::vect), orientation groups (@f$\SO n@f$, see SOn.h).
Furthermore it provides a convenient way to define compound manifolds such as poses
@f$ \SE n \simeq \SO n \times \R^n@f$), with components being accessible by name instead of index. 

All manifolds are accessible by two methods
\code
  void boxplus(MTK::vectview<const double, DOF> delta, double scale);
  void boxminus(MTK::vectview<double, DOF> result, const Manifold &oth);
\endcode
where for a @c Manifold m @c boxplus stores @f$m \mplus scale\cdot delta@f$ to m
and @c boxminus stores @f$m \mmnus oth \to result @f$.


For matrix/vector operations, MTK uses the open source library
Eigen 2 (http://eigen.tuxfamily.org) making the implementation of measurement
functions straight-forward using a MATLAB-like notation
via overloaded operators, and efficient using expression templates.

*/


/**
\page MANIFOLD Introduction to Manifolds

A manifold is a structure which locally looks like a Euclidean vector space (@f$\R^n@f$)
but in general does not globally. For a first more detailed introduction see the 
<a href="http://en.wikipedia.org/Manifold"> wikipedia article on manifolds </a> 
or read J.M. Lee's book "Introduction to Smooth Manifolds" \cite{Lee2003}.

In MTK a manifold @f$\S@f$ is encapsulated using two operators 
@f{alignat*}{2}
\mplus &: \S\times\R^n &&\to \S, \\
\mmnus &: \S\times \S  &&\to\R^n. 
@f}
The @f$\mplus@f$ operator adds (preferably small) perturbation to a state, whereas the 
@f$\mmnus@f$ operator is its inverse, i.e. it calculates the perturbation necessary to
come from on state to another:
@f[
x \mplus(y\mmnus x) = y.
@f]
The operators are implemented as memberfunctions @c boxplus and @c boxminus as shown in @ref MTKTutorial.
For more technical details see our paper \cite{Hertzberg2011}.
*/


/**
\page BUILD Building MTK and SLoM

MTK is completely template-based, so no building is required. 
Building projects using MTK requires the MTK as well as Eigen (version 2.0.9 or later -- 
Eigen 3 is not supported at the moment) and the Boost preprocessor library in the include path.
*/
