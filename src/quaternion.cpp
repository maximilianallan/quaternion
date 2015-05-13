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

Quaternion::Quaternion() : internal_quaternion_(1, 0, 0, 0) {

}

Quaternion::Quaternion(const boost::math::quaternion<double> &x) :internal_quaternion_(x) {

}


std::ostream &operator<<(std::ostream &stream, const Quaternion &a){
  stream << a.W() << ", " << a.X() << ", " << a.Y() << ", " << a.Z();
  return stream;
}

Quaternion::Quaternion(const cv::Vec3d &euler_angles){

  const double half_phi = euler_angles[0] / 2;
  const double half_theta = euler_angles[1] / 2;
  const double half_psi = euler_angles[2] / 2;
  const double q1 = (cos(half_phi)*cos(half_theta)*cos(half_psi)) + (sin(half_phi)*sin(half_theta)*sin(half_psi));
  const double q2 = (sin(half_phi)*cos(half_theta)*cos(half_psi)) - (cos(half_phi)*sin(half_theta)*sin(half_psi));
  const double q3 = (cos(half_phi)*sin(half_theta)*cos(half_psi)) + (sin(half_phi)*cos(half_theta)*sin(half_psi));
  const double q4 = (cos(half_phi)*cos(half_theta)*sin(half_psi)) - (sin(half_phi)*sin(half_theta)*cos(half_psi));

  internal_quaternion_ = boost::math::quaternion<double>(q1, q2, q3, q4);
}

Quaternion::Quaternion(const cv::Mat &rotation_matrix) {
  if (rotation_matrix.size() != cv::Size(3, 3) && rotation_matrix.type() != CV_64FC1) throw(std::runtime_error("Error, rotation matrix used to construct quaternion should be 3x3 double matrix!\n"));
  double w = sqrt(1.0 + rotation_matrix.at<double>(0, 0) + rotation_matrix.at<double>(1, 1) + rotation_matrix.at<double>(2, 2)) / 2.0;
  double w4 = 4.0*w;
  double x = (rotation_matrix.at<double>(2, 1) - rotation_matrix.at<double>(1, 2)) / w4;
  double y = (rotation_matrix.at<double>(0, 2) - rotation_matrix.at<double>(2, 0)) / w4;
  double z = (rotation_matrix.at<double>(1, 0) - rotation_matrix.at<double>(0, 1)) / w4;
  internal_quaternion_ = boost::math::quaternion<double>(w, x, y, z);

}

cv::Vec3d Quaternion::EulerAngles() const {
  const double phi = (double)atan2(2 * (W()*X() + Y()*Z()), 1 - 2 * (X()*X() + Y()*Y()));
  const double theta = (double)asin(2 * (W()*Y() - Z()*X()));
  const double psi = (double)atan2(2 * (W()*Z() + X()*Y()), 1 - 2 * (Y()*Y() + Z()*Z()));
  return cv::Vec3d(phi, theta, psi);
  //return cv::Vec3d(theta,psi,phi)
}

cv::Vec3d Quaternion::AngleAxis() const {

  //const double theta = 2 * acos( W() );
  const double l2_norm = (double)sqrt((X()*X()) + (Y()*Y()) + (Z()*Z()));
  const double theta = (double)(2 * atan2((double)l2_norm, W()));
  if (theta == 0.0f) return cv::Vec3d(0, 0, 0);
  const cv::Vec3d omega((double)X() / l2_norm, (double)Y() / l2_norm, (double)Z() / l2_norm);

  return theta * omega;

}

cv::Mat Quaternion::RotationMatrix() const {

  cv::Mat ret(3, 3, CV_64FC1);

  ret.at<double>(0, 0) = (1 - (2 * Y()*Y()) - (2 * Z()*Z()));
  ret.at<double>(0, 1) = (2 * X()*Y()) - (2 * Z()*W());
  ret.at<double>(0, 2) = (2 * X()*Z()) + (2 * Y()*W());

  ret.at<double>(1, 0) = (2 * X()*Y()) + (2 * Z()*W());
  ret.at<double>(1, 1) = 1 - (2 * X()*X()) - (2 * Z()*Z());
  ret.at<double>(1, 2) = (2 * Y()*Z()) - (2 * X()*W());

  ret.at<double>(2, 0) = (2 * X()*Z()) - (2 * Y()*W());
  ret.at<double>(2, 1) = (2 * Y()*Z()) + (2 * X()*W());
  ret.at<double>(2, 2) = 1 - (2 * X()*X()) - (2 * Y()*Y());

  return ret;

}


Quaternion Quaternion::Inverse() const {
  return Quaternion(boost::math::conj<double>(internal_quaternion_));
}
double Quaternion::AngularDistanceToQuaternion(const Quaternion &other){

  /* FLOATING POINT ERRORS - this is not feasible
  if ( this->Normalize() != *this ) {
  std::cerr << *this << "\n" << this->Normalize() << "\n";
  throw(std::runtime_error("Error, quaternion caller is not normalized!\n"));
  }
  if ( other.Normalize() != other ) {
  std::cerr << other << "\n" << other.Normalize() << "\n";
  throw(std::runtime_error("Error, other quaternion is not normalized!\n"));
  }*/

  double inner_sq = X()*other.X() + Y()*other.Y() + Z()*other.Z() + W()*other.W();
  inner_sq = inner_sq * inner_sq;
  return acos(2 * (inner_sq)-1);
}

Quaternion::Quaternion(const double angle, const cv::Vec3d &axis){
  double norm = 0;
  for (int i = 0; i<3; i++) norm += axis[i] * axis[i];
  norm = std::sqrt(norm);
  if (norm == 0) norm = 0.00001;
  cv::Vec3d axis_normed;
  for (int i = 0; i<3; i++) axis_normed = axis[i] / norm;
  const double sin_angle_2 = sin(angle / 2);
  internal_quaternion_ = boost::math::quaternion<double>(cos(angle / 2), axis[0] * sin_angle_2, axis[1] * sin_angle_2, axis[2] * sin_angle_2);
  *this = Normalize();
}

Quaternion Quaternion::FromVectorToVector(const cv::Vec3d &from, const cv::Vec3d to){

  cv::Vec3d from_n, to_n;
  cv::normalize(from, from_n);
  cv::normalize(to, to_n);

  double d = from_n.dot(to_n);
  if (d >= 1.0){
    //return identity quaternion
    return boost::math::quaternion<double>(1, 0, 0, 0);
  }

  if (d <= -1.0){
    //for cases when the vectors point in opposite directions
    //see http://irrlicht.sourceforge.net/docu/quaternion_8h_source.html line 658
    cv::Vec3d ax(1, 0, 0);
    ax = ax.cross(from_n);
    if (ax.dot(ax) == 0){
      ax = cv::Vec3d(0, 1, 0);
      ax = ax.cross(from_n);
    }
    Quaternion qr(0, ax);
    return qr.Normalize();
  }

  double s = (double)sqrt((1 + d) * 2);
  double inv_s = 1.0f / s;

  cv::Vec3d axis = from_n.cross(to_n);
  Quaternion q(boost::math::quaternion<double>(s*0.5f, axis[0] * inv_s, axis[1] * inv_s, axis[2] * inv_s));
  return q.Normalize();

}

Quaternion Quaternion::Normalize() const {

  //Note boost::math::norm(quaternion) is Cayley norm NOT Euclidean norm...
  double mag_2 = (internal_quaternion_.R_component_1() * internal_quaternion_.R_component_1())
    + (internal_quaternion_.R_component_2() * internal_quaternion_.R_component_2())
    + (internal_quaternion_.R_component_3() * internal_quaternion_.R_component_3())
    + (internal_quaternion_.R_component_4() * internal_quaternion_.R_component_4());

  double mag = sqrt(mag_2);

  return boost::math::quaternion<double>(internal_quaternion_.R_component_1() / mag,
    internal_quaternion_.R_component_2() / mag,
    internal_quaternion_.R_component_3() / mag,
    internal_quaternion_.R_component_4() / mag);
}

cv::Vec3d Quaternion::RotateVector(const cv::Vec3d &to_rotate) const {

  cv::Vec3d v = to_rotate;
  cv::Vec3d uv, uuv;
  cv::Vec3d qvec((double)this->X(), (double)this->Y(), (double)this->Z());
  uv = qvec.cross(v);
  uuv = qvec.cross(uv);
  uv *= (2.0f * this->W());
  uuv *= 2.0f;

  return v + uv + uuv;

  /*
  const boost::math::quaternion<double> vec_quat(0,to_rotate[0],to_rotate[1],to_rotate[2]);

  boost::math::quaternion<double> rotated = (internal_quaternion_ * vec_quat) * boost::math::conj<double>(internal_quaternion_);

  return cv::Vec3d((double)rotated.R_component_2(),(double)rotated.R_component_3(),(double)rotated.R_component_4());
  */
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