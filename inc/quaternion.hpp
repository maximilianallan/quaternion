/**
Copyright (c) 2013, MAXIMILIAN ALLAN
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met: 

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies, 
either expressed or implied, of the FreeBSD Project.
*/

#ifndef __QUATERNION_HPP__
#define __QUATERNION_HPP__

#include <cv.h>
#include <boost/math/quaternion.hpp>

namespace sv {

 class Quaternion {
  
  public:

    Quaternion();
    Quaternion(const boost::math::quaternion<double> &x);
    Quaternion(const double angle, const cv::Vec3d &axis);
    explicit Quaternion(const cv::Vec3d &euler_angles);
    explicit Quaternion(const cv::Mat &rotation_matrix);
    static Quaternion FromVectorToVector(const cv::Vec3d &from, const cv::Vec3d to);
    
    
    Quaternion Normalize() const ;
    inline Quaternion Inverse() const;
    inline double X() const;
    inline double Y() const;
    inline double Z() const;
    inline double W() const;

    //(psi,theta,phi)
    inline cv::Mat RotationMatrix() const ;
    inline cv::Vec3d EulerAngles() const ;
    inline cv::Vec3d AngleAxis() const ;
    inline cv::Vec3d RotateVector(const cv::Vec3d &to_rotate) const ;
    inline double AngularDistanceToQuaternion(const Quaternion &other);
    
    inline friend std::ostream &operator<<(std::ostream &stream, const Quaternion &a);
    friend Quaternion operator*(const Quaternion &a, const Quaternion &b);
    friend Quaternion operator+(const Quaternion &a, const Quaternion &b);
    friend Quaternion operator-(const Quaternion &a, const Quaternion &b);
    friend bool operator==(const Quaternion &a, const Quaternion &b);
    friend bool operator!=(const Quaternion &a, const Quaternion &b);

  protected:

    boost::math::quaternion<double> internal_quaternion_;


  };



 Quaternion operator-(const Quaternion &a, const Quaternion &b) {

   Quaternion q;
   q.internal_quaternion_ = a.internal_quaternion_ - b.internal_quaternion_;
   return q;

 }

 bool operator==(const Quaternion &a, const Quaternion &b){

   return a.X() == b.X() && a.Y() == b.Y() && a.Z() == b.Z() && a.W() == b.W();

 }

 bool operator!=(const Quaternion &a, const Quaternion &b){
   return !(a == b);
 }


 Quaternion operator+(const Quaternion &a, const Quaternion &b) {

   Quaternion q;
   q.internal_quaternion_ = a.internal_quaternion_ + b.internal_quaternion_;
   return q;

 }


 Quaternion operator*(const Quaternion &a, const Quaternion &b) {

   Quaternion q;
   q.internal_quaternion_ = a.internal_quaternion_ * b.internal_quaternion_;
   return q;
 }






  
}

#endif
