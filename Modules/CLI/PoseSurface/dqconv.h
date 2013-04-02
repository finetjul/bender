/* dqconv.h

  Conversion routines between (regular quaternion, translation) and dual quaternion.

  Version 1.0.0, February 7th, 2007

  Copyright (C) 2006-2007 University of Dublin, Trinity College, All Rights
  Reserved

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the author(s) be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  Author: Ladislav Kavan, kavanl@cs.tcd.ie
  Changes:
    - changed float to double for convenience (Yuanxin Liu)

*/

#ifndef __dqconv_h
#define __dqconv_h

#include <math.h>

#include <vtkDualQuaternion.h>

// input: unit quaternion 'q0', translation vector 't'
// output: unit dual quaternion 'dq'
void QuatTrans2UDQ(const double q0[4], const double t[3],
                  double dq[2][4])
{
   // non-dual part (just copy q0):
   for (int i=0; i<4; i++) dq[0][i] = q0[i];
   // dual part:
   dq[1][0] = -0.5*(t[0]*q0[1] + t[1]*q0[2] + t[2]*q0[3]);
   dq[1][1] = 0.5*( t[0]*q0[0] + t[1]*q0[3] - t[2]*q0[2]);
   dq[1][2] = 0.5*(-t[0]*q0[3] + t[1]*q0[0] + t[2]*q0[1]);
   dq[1][3] = 0.5*( t[0]*q0[2] - t[1]*q0[1] + t[2]*q0[0]);
}

// input: unit dual quaternion 'dq'
// output: unit quaternion 'q0', translation vector 't'
void UDQ2QuatTrans(const double dq[2][4],
                  double q0[4], double t[3])
{
   // regular quaternion (just copy the non-dual part):
   for (int i=0; i<4; i++) q0[i] = dq[0][i];
   // translation vector:
   t[0] = 2.0*(-dq[1][0]*dq[0][1] + dq[1][1]*dq[0][0] - dq[1][2]*dq[0][3] + dq[1][3]*dq[0][2]);
   t[1] = 2.0*(-dq[1][0]*dq[0][2] + dq[1][1]*dq[0][3] + dq[1][2]*dq[0][0] - dq[1][3]*dq[0][1]);
   t[2] = 2.0*(-dq[1][0]*dq[0][3] - dq[1][1]*dq[0][2] + dq[1][2]*dq[0][1] + dq[1][3]*dq[0][0]);
}

// input: dual quat. 'dq' with non-zero non-dual part
// output: unit quaternion 'q0', translation vector 't'
void DQ2QuatTrans(const double dq[2][4],
                  double q0[4], double t[3])
{
   double len = 0.0;
   for (int i=0; i<4; i++) len += dq[0][i] * dq[0][i];
   len = sqrt(len);
   for (int i=0; i<4; i++) q0[i] = dq[0][i] / len;
   t[0] = 2.0*(-dq[1][0]*dq[0][1] + dq[1][1]*dq[0][0] - dq[1][2]*dq[0][3] + dq[1][3]*dq[0][2]) / len;
   t[1] = 2.0*(-dq[1][0]*dq[0][2] + dq[1][1]*dq[0][3] + dq[1][2]*dq[0][0] - dq[1][3]*dq[0][1]) / len;
   t[2] = 2.0*(-dq[1][0]*dq[0][3] - dq[1][1]*dq[0][2] + dq[1][2]*dq[0][1] + dq[1][3]*dq[0][0]) / len;
}

void ScLERP(const double a[2][4], const double b[2][4], double t, double out[2][4])
{
  vtkDualQuaternion<double> q1(a[0][0], a[0][1], a[0][2], a[0][3],
                               a[1][0], a[1][1], a[1][2], a[1][3]);
  vtkDualQuaternion<double> q2(b[0][0], b[0][1], b[0][2], b[0][3],
                               b[1][0], b[1][1], b[1][2], b[1][3]);
  vtkDualQuaternion<double> q3 = q1.ScLerp2(t, q2);
  out[0][0] = q3.GetReal().GetW();
  out[0][1] = q3.GetReal().GetX();
  out[0][2] = q3.GetReal().GetY();
  out[0][3] = q3.GetReal().GetZ();
  out[1][0] = q3.GetDual().GetW();
  out[1][1] = q3.GetDual().GetX();
  out[1][2] = q3.GetDual().GetY();
  out[1][3] = q3.GetDual().GetZ();
}

#endif
