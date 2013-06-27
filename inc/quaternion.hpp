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

    Quaternion(){}
    Quaternion(const boost::math::quaternion<double> &x):internal_quaternion_(x) {}
    Quaternion(const double angle, const cv::Vec3f &axis);
    inline explicit Quaternion(const cv::Vec3f &euler_angles);
    static Quaternion FromVectorToVector(const cv::Vec3f &from, const cv::Vec3f to);
    cv::Vec3f RotateVector(const cv::Vec3f &to_rotate) const ;
    Quaternion Normalize() const ;
    inline Quaternion Inverse() const;
    inline double X() const;
    inline double Y() const;
    inline double Z() const;
    inline double W() const;

    //(psi,theta,phi)
    inline cv::Vec3f EulerAngles() const ;
    inline friend std::ostream &operator<<(std::ostream &stream, const Quaternion &a);
    inline friend Quaternion operator*(const Quaternion &a, const Quaternion &b);
    inline friend Quaternion operator+(const Quaternion &a, const Quaternion &b);
    inline friend Quaternion operator-(const Quaternion &a, const Quaternion &b);

  protected:

    boost::math::quaternion<double> internal_quaternion_;


  };

  std::ostream &operator<<(std::ostream &stream, const Quaternion &a){
    stream << "[" << a.X() << ", " << a.X() << ", " << a.Y() << ", " << a.Z() << "]\n";
    return stream;
  }

  Quaternion::Quaternion(const cv::Vec3f &euler_angles){

    const float half_phi = euler_angles[0]/2;
    const float half_theta = euler_angles[1]/2;
    const float half_psi = euler_angles[2]/2;
    const float q1 = (cos(half_phi)*cos(half_theta)*cos(half_psi)) + (sin(half_phi)*sin(half_theta)*sin(half_psi));
    const float q2 = (sin(half_phi)*cos(half_theta)*cos(half_psi)) - (cos(half_phi)*sin(half_theta)*sin(half_psi));
    const float q3 = (cos(half_phi)*sin(half_theta)*cos(half_psi)) + (sin(half_phi)*cos(half_theta)*sin(half_psi));
    const float q4 = (cos(half_phi)*cos(half_theta)*sin(half_psi)) - (sin(half_phi)*sin(half_theta)*cos(half_psi));
    
    internal_quaternion_ = boost::math::quaternion<double>(q1,q2,q3,q4);
  }

  cv::Vec3f Quaternion::EulerAngles() const {
    const float phi = atan2( 2*(W()*X() + Y()*Z()), 1-2*(X()*X()+Y()*Y()) );
    const float theta = asin(2*(W()*Y() - Z()*X()));
    const float psi = atan2( 2*(W()*Z() + X()*Y()), 1-2*(Y()*Y()+Z()*Z()) );
    return cv::Vec3f(phi,theta,psi);
  }

  Quaternion Quaternion::Inverse() const {
    return Quaternion(boost::math::conj<double>(internal_quaternion_));
  }
  
  
  Quaternion operator-(const Quaternion &a, const Quaternion &b) {

    Quaternion q;
    q.internal_quaternion_ = a.internal_quaternion_ - b.internal_quaternion_;
    return q;
  
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

  double Quaternion::X() const {
    return internal_quaternion_.R_component_2();
  }
  double Quaternion::Y() const {
    return internal_quaternion_.R_component_3();
  }
  double Quaternion::Z() const {
    return internal_quaternion_.R_component_4();
  }
  double Quaternion::W() const {
    return internal_quaternion_.R_component_1();
  }
}

#endif
