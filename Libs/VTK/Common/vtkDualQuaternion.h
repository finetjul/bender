/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkQuaternion.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef __vtkDualQuaternion_h
#define __vtkDualQuaternion_h

#include "vtkQuaternion.h"

template<typename T> class vtkDualQuaternion
{
public:
  vtkDualQuaternion();
  vtkDualQuaternion(const T& rw, const T& rx, const T& ry, const T& rz,
                    const T& dw, const T& dx, const T& dy, const T& dz);
  vtkDualQuaternion(const vtkQuaternion<T>& real, const vtkQuaternion<T>& dual);

  vtkQuaternion<T> GetReal()const;
  vtkQuaternion<T> GetDual()const;

  ///
  void Invert();
  /// Conjugate / SquaredNorm
  vtkDualQuaternion Inverse()const;
  vtkDualQuaternion Inverse2()const;
  void Normalize();
  void Conjugate();
  vtkDualQuaternion<T> Conjugated()const;

  void SetRotationPosition(const vtkQuaternion<T>& rotation, const T position[3]);
  void SetRotation(const vtkQuaternion<T>& rotation);
  void SetTranslation(const T translation[3]);
  void GetPosition(T position[3])const;

  bool operator==(const vtkDualQuaternion<T>& dq)const;
  vtkDualQuaternion<T> operator+(const vtkDualQuaternion<T>& dq)const;
  vtkDualQuaternion<T> operator-(const vtkDualQuaternion<T>& dq)const;
  vtkDualQuaternion<T> operator-()const;
  vtkDualQuaternion<T> operator*(const T& scalar)const;
  vtkDualQuaternion<T> operator*(const vtkDualQuaternion<T>& dq)const;

  vtkDualQuaternion<T> Lerp(const T& t, const vtkDualQuaternion<T>& dq)const;
  vtkDualQuaternion<T> ScLerp(const T& t, const vtkDualQuaternion<T>& dq)const;
  vtkDualQuaternion<T> ScLerp2(const T& t, const vtkDualQuaternion<T>& dq)const;
  vtkDualQuaternion<T> Dot(const vtkDualQuaternion<T>& dq)const;
  void LengthSquared(T& real, T& dual)const;
  void ReciprocalLengthSquared(T& real, T& dual)const;
  void ToScrew(T& angle, T& pitch, T dir[3], T moment[3])const;
  void FromScrew(T angle, T pitch, T dir[3], T moment[3]);

protected:
  vtkQuaternion<T> Real;
  vtkQuaternion<T> Dual;
};

#include "vtkDualQuaternion.txx"

#endif
