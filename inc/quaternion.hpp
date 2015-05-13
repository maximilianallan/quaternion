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

#include <cmath>
#include <iostream>
#include <cv.h>

namespace math {

  const static double eps = 1e-9;
  
  class Quaternion {
    
  public:
    
    Quaternion();
    Quaternion(const double w, const double x, const double y, const double z);
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
    
    friend inline std::ostream &operator<<(std::ostream &stream, const Quaternion &a);
    Quaternion operator*(const Quaternion &rhs);
    Quaternion operator+(const Quaternion &rhs);
    Quaternion operator-(const Quaternion &rhs);
    bool operator==(const Quaternion &rhs);
    bool operator!=(const Quaternion &rhs);
    
  protected:

    template<typename T>
    void InitFromMatrix(const cv::Mat &matrix);
    
    double w_;
    double x_;
    double y_;
    double z_;
    
  };
  
  template<typename T>
  void Quaternion::InitFromMatrix(const cv::Mat &rotation_matrix){
    
    w_ = sqrt(1.0 + rotation_matrix.at<T>(0, 0) + rotation_matrix.at<T>(1, 1) + rotation_matrix.at<T>(2, 2)) / 2.0;
    const double w4 = 4.0*w_;
    x_ = (rotation_matrix.at<T>(2, 1) - rotation_matrix.at<T>(1, 2)) / w4;
    y_ = (rotation_matrix.at<T>(0, 2) - rotation_matrix.at<T>(2, 0)) / w4;
    z_ = (rotation_matrix.at<T>(1, 0) - rotation_matrix.at<T>(0, 1)) / w4;

  }
 
}

#endif
