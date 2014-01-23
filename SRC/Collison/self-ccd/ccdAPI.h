/*************************************************************************\

  Copyright 2010 The University of North Carolina at Chapel Hill.
  All Rights Reserved.

  Permission to use, copy, modify and distribute this software and its
  documentation for educational, research and non-profit purposes, without
   fee, and without a written agreement is hereby granted, provided that the
  above copyright notice and the following three paragraphs appear in all
  copies.

  IN NO EVENT SHALL THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL BE
  LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
  CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE
  USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY
  OF NORTH CAROLINA HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH
  DAMAGES.

  THE UNIVERSITY OF NORTH CAROLINA SPECIFICALLY DISCLAIM ANY
  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE
  PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
  NORTH CAROLINA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,
  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

  The authors may be contacted via:

  US Mail:             GAMMA Research Group at UNC
                       Department of Computer Science
                       Sitterson Hall, CB #3175
                       University of N. Carolina
                       Chapel Hill, NC 27599-3175

  Phone:               (919)962-1749

  EMail:              geom@cs.unc.edu; tang_m@zju.edu.cn


\**************************************************************************/

#ifndef _CCDAPI_H_
#define _CCDAPI_H_

#include "vec3f.h"
#include "feature.h"

  typedef void ccdEETestCallback(unsigned int e1_v1, unsigned int e1_v2, unsigned int e2_v1, unsigned int e2_v2, float t);

/** 
 * A function which will be called when a face-vertex collision is detected.
 * 
 * @param vid the index of the vertex in collision.
 * @param fid the index of the face in collision.
 * @param t a time in [0,1] indicates when the collision happens.
 */
typedef void ccdVFTestCallback(unsigned int vid, unsigned int fid, float t);

extern void ccdInitModel(SELF_CCD::vec3f_list &, SELF_CCD::tri_list &);
extern void ccdUpdateVtxs(SELF_CCD::vec3f_list &);
extern void ccdQuitModel();
extern void ccdChecking(bool);
extern void ccdReport();
extern void ccdSetEECallback(ccdEETestCallback *funcEE);
extern void ccdSetVFCallback(ccdVFTestCallback *funcVF);

#endif /* _CCDAPI_H_ */
