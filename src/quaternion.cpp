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

#include "../inc/quaternion.hpp"


using namespace sv;

Quaternion::Quaternion(const double angle, const cv::Vec3d &axis):internal_quaternion_(angle,axis[0],axis[1],axis[2]){}

Quaternion Quaternion::FromVectorToVector(const cv::Vec3d &from, const cv::Vec3d to){
  
  cv::Vec3d from_n,to_n;
  cv::normalize(from,from_n);
  cv::normalize(to,to_n);
      
  double d = from_n.dot(to_n);

  if(d >= 1.0){
    return boost::math::quaternion<double>();
  }

  //check if d \approx = 0

  double s = (double)sqrt( (1+d)*2 );
  double inv_s = 1.0f/s;

  cv::Vec3d axis = from_n.cross(to_n);

  Quaternion q( s*0.5f, cv::Vec3d(axis[0]*inv_s, axis[1]*inv_s, axis[2]*inv_s ));

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

cv::Vec3d Quaternion::RotateVector(const cv::Vec3d &to_rotate) const {
   
  const boost::math::quaternion<double> vec_quat(0,to_rotate[0],to_rotate[1],to_rotate[2]);

  boost::math::quaternion<double> rotated = (internal_quaternion_ * vec_quat) * boost::math::conj<double>(internal_quaternion_);

  return cv::Vec3d((double)rotated.R_component_2(),(double)rotated.R_component_3(),(double)rotated.R_component_3());

}
