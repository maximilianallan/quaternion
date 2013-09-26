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
    inline Quaternion(const double angle, const cv::Vec3f &axis);
    inline explicit Quaternion(const cv::Vec3f &euler_angles);
    inline explicit Quaternion(const cv::Mat &rotation_matrix);
    inline static Quaternion FromVectorToVector(const cv::Vec3f &from, const cv::Vec3f to);
    inline cv::Vec3f RotateVector(const cv::Vec3f &to_rotate) const ;
    inline Quaternion Normalize() const ;
    inline Quaternion Inverse() const;
    inline double X() const;
    inline double Y() const;
    inline double Z() const;
    inline double W() const;

    //(psi,theta,phi)
    inline cv::Vec3f EulerAngles() const ;
    inline cv::Vec3f AngleAxis() const ;
    inline friend std::ostream &operator<<(std::ostream &stream, const Quaternion &a);
    inline friend Quaternion operator*(const Quaternion &a, const Quaternion &b);
    inline friend Quaternion operator+(const Quaternion &a, const Quaternion &b);
    inline friend Quaternion operator-(const Quaternion &a, const Quaternion &b);

  protected:

    boost::math::quaternion<double> internal_quaternion_;


  };

  std::ostream &operator<<(std::ostream &stream, const Quaternion &a){
    stream << a.W() << " " << a.X() << " " << a.Y() << " " << a.Z();
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

  Quaternion::Quaternion(const cv::Mat &rotation_matrix) {
    if (rotation_matrix.size() != cv::Size(3,3) && rotation_matrix.type() != CV_64FC1) throw(std::runtime_error("Error, rotation matrix used to construct quaternion should be 3x3 double matrix!\n"));
    double w = sqrt( 1.0 + rotation_matrix.at<double>(0,0) + rotation_matrix.at<double>(1,1) + rotation_matrix.at<double>(2,2))/2.0;
    double w4 = 4.0*w;
    double x = (rotation_matrix.at<double>(2,1) - rotation_matrix.at<double>(1,2))/w4;
    double y = (rotation_matrix.at<double>(0,2) - rotation_matrix.at<double>(2,0))/w4;
    double z = (rotation_matrix.at<double>(1,0) - rotation_matrix.at<double>(0,1))/w4;
    internal_quaternion_ = boost::math::quaternion<double>(w,x,y,z);

  }

  cv::Vec3f Quaternion::EulerAngles() const {
    const float phi = (float)atan2( 2*(W()*X() + Y()*Z()), 1-2*(X()*X()+Y()*Y()) );
    const float theta = (float)asin(2*(W()*Y() - Z()*X()));
    const float psi = (float)atan2( 2*(W()*Z() + X()*Y()), 1-2*(Y()*Y()+Z()*Z()) );
    return cv::Vec3f(phi,theta,psi);
    //return cv::Vec3f(theta,psi,phi)
  }

  cv::Vec3f Quaternion::AngleAxis() const {

    //const float theta = 2 * acos( W() );
    const float l2_norm = sqrt( (X()*X()) + (Y()*Y()) + (Z()*Z()) );
    const float theta = 2 * atan2( (double)l2_norm , W() );
    if(theta == 0.0f) return cv::Vec3f(0,0,0);
    const cv::Vec3f omega( X()/l2_norm, Y()/l2_norm, Z()/l2_norm );

    return omega;

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
  
  
Quaternion::Quaternion(const double angle, const cv::Vec3f &axis):internal_quaternion_(angle,axis[0],axis[1],axis[2]){}

Quaternion Quaternion::FromVectorToVector(const cv::Vec3f &from, const cv::Vec3f to){
  
  cv::Vec3f from_n,to_n;
  cv::normalize(from,from_n);
  cv::normalize(to,to_n);
    
  float d = from_n.dot(to_n);
  if(d >= 1.0){
    //return identity quaternion
    return boost::math::quaternion<double>(1,0,0,0);
  }

  if(d <= -1.0){
    //for cases when the vectors point in opposite directions
    //see http://irrlicht.sourceforge.net/docu/quaternion_8h_source.html line 658
    cv::Vec3f ax(1,0,0);
    ax = ax.cross(from_n);
    if(ax.dot(ax) == 0){
      ax = cv::Vec3f(0,1,0);
      ax = ax.cross(from_n);
    }
    Quaternion qr(0,ax);
    return qr.Normalize();
  }

  float s = (float)sqrt( (1+d)*2 );
  float inv_s = 1.0f/s;

  cv::Vec3f axis = from_n.cross(to_n);
  Quaternion q( s*0.5f, cv::Vec3f(axis[0]*inv_s, axis[1]*inv_s, axis[2]*inv_s ));
  return q.Normalize();

}

Quaternion Quaternion::Normalize() const {
  
  //Note boost::math::norm(quaternion) is Cayley norm NOT Euclidean norm...
  double mag_2 = (internal_quaternion_.R_component_1() * internal_quaternion_.R_component_1()) 
    + (internal_quaternion_.R_component_2() * internal_quaternion_.R_component_2()) 
    + (internal_quaternion_.R_component_3() * internal_quaternion_.R_component_3()) 
    + (internal_quaternion_.R_component_4() * internal_quaternion_.R_component_4());
      
  double mag = sqrt(mag_2);
      
  return boost::math::quaternion<double>(internal_quaternion_.R_component_1()/mag,
                                         internal_quaternion_.R_component_2()/mag,
                                         internal_quaternion_.R_component_3()/mag,
                                         internal_quaternion_.R_component_4()/mag);
}

cv::Vec3f Quaternion::RotateVector(const cv::Vec3f &to_rotate) const {
   
  const boost::math::quaternion<double> vec_quat(0,to_rotate[0],to_rotate[1],to_rotate[2]);

  boost::math::quaternion<double> rotated = (internal_quaternion_ * vec_quat) * boost::math::conj<double>(internal_quaternion_);

  return cv::Vec3f((float)rotated.R_component_2(),(float)rotated.R_component_3(),(float)rotated.R_component_4());

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
