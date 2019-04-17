#include <iostream>
#include <vector>
#include <Eigen/Dense>

typedef std::pair<Eigen::Vector3d, double> returnType; //shorter

Eigen::Matrix3d Euler2A(double fi, double teta, double psi){
    Eigen::Matrix3d rx, ry, rz;
    rx << 1.0, 0.0,           0.0,
          0.0, std::cos(fi), -std::sin(fi),
          0.0, std::sin(fi),  std::cos(fi);

    ry <<  std::cos(teta), 0.0, std::sin(teta),
           0.0,            1.0, 0.0,
          -std::sin(teta), 0.0, std::cos(teta);

    rz << std::cos(psi), -std::sin(psi), 0.0,
          std::sin(psi),  std::cos(psi), 0.0,
          0.0,            0.0,           1.0;

 return rz*ry*rx;
}

returnType axisAngle(const Eigen::Matrix3d &A){

    Eigen::Matrix3d E = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d tmp = A - E;
    Eigen::Vector3d p = tmp.row(0).cross(tmp.row(1));
    p.normalize();

    Eigen::Vector3d u = tmp.row(1);
    Eigen::Vector3d up = A*u;

    double teta = std::acos(u.dot(up)/(u.norm()*up.norm()));

    Eigen::Matrix3d mtp;
    mtp << u, up, p;
    if(mtp.determinant() < 0)
        p *= -1;

 return returnType(p.transpose(),teta);
}

Eigen::Matrix3d Rodrigez(Eigen::Vector3d p, double teta){

    Eigen::Matrix3d px;
    px <<  0.0,  -p(2),  p(1),
           p(2),  0.0,  -p(0),
          -p(1),  p(0),  0.0;

    Eigen::Matrix3d E = Eigen::Matrix3d::Identity();

 return p*p.transpose() + std::cos(teta)*(E-p*p.transpose()) + std::sin(teta)*px;
}

Eigen::Vector3d A2Euler(Eigen::Matrix3d A){

    double fi, teta, psi;

    if(A(2,0) < 1){
        if(A(2,0) > -1){
            psi = std::atan2(A(1,0), A(0,0));
            teta = std::asin(-A(2,0));
            fi = std::atan2(A(2,1), A(2,2));
        }
        else{
            psi = std::atan2(-A(0,1), A(1,1));
            teta = M_PI/2;
            fi = 0.0;
        }
    }
    else{
        psi = std::atan2(-A(0,1), A(1,1));
        teta = -M_PI/2;
        fi = 0.0;
    }

 return Eigen::Vector3d{fi, teta, psi};
}

Eigen::Vector4d AxisAngle2Q(Eigen::Vector3d p, double teta){

    double w = std::cos(teta/2);
    p.normalize();
    Eigen::Vector3d xyz = std::sin(teta/2) * p;

 return Eigen::Vector4d{xyz(0), xyz(1), xyz(2), w};
}

returnType Q2AxisAngle(Eigen::Vector4d q){

    q.normalize();
    double teta = 2*std::acos(q(3)); //2*arccos(w)
    Eigen::Vector3d v;
    Eigen::Vector3d tmp{q(0), q(1), q(2)};
    if(std::abs(q(3)) == 1.0)
        v << 1.0, 0.0, 0.0;
    else
        v = tmp.normalized();

 return returnType(v.transpose(), teta);
}

Eigen::Vector4d Slerp(Eigen::Vector4d q1, Eigen::Vector4d q2, double time, double currentTime) {
    double cos = q1.transpose()*q2;
    if(cos < 0) {
        q1 *= -1;
        cos *= -1;
    }

    if(cos > 0.95) {
        return q1;
    }

    double fi = std::acos(cos);
    double first = std::sin(fi * (1 - currentTime / time)) / std::sin(fi);
    double second = std::sin(fi * currentTime / time) / std::sin(fi);

 return first*q1 + second*q2;
}

int main() {

    Eigen::Matrix3d A = Euler2A(std::asin(2.0/3), -std::acos(2.0/5), std::atan(1.0/2));
    std::cout << "Euler2A[φ, θ, ψ]: \n";
    std::cout << "A:\n" << A << "\n";
    std::cout << "---------------------------------------------------" << std::endl;

    returnType rt = axisAngle(A);
    Eigen::Vector3d p = rt.first;
    double teta = rt.second;
    std::cout << "AxisAngle[A]\n";
    std::cout << "p:\n" << p << std::endl << "teta:\n " << teta << "\n";
    std::cout << "---------------------------------------------------" << std::endl;

    Eigen::Matrix3d Ap = Rodrigez(p, teta);
    std::cout << "Rodrigez[p, φ]\n";
    std::cout << "A:\n" << Ap << "\n";
    std::cout << "---------------------------------------------------" << std::endl;

    Eigen::Vector3d angles = A2Euler(A);
    std::cout << "A2Euler[A]\n";
    std::cout << "Euler's angles:\n " << angles.transpose() << "\n";
    std::cout << "---------------------------------------------------" << std::endl;

    Eigen::Vector4d quat = AxisAngle2Q(p, teta);
    std::cout << "AxisAngle2Q[p, φ]\n";
    std::cout << "Quat:\n" << quat.transpose() << "\n";
    std::cout << "---------------------------------------------------" << std::endl;

    returnType qrt = Q2AxisAngle(quat);
    Eigen::Vector3d qp = qrt.first;
    double qteta = qrt.second;
    std::cout << "Q2AxisAngle[q]\n";
    std::cout << "p:\n" << qp << std::endl << "teta:\n " << qteta << "\n" ;

 return 0;
}
