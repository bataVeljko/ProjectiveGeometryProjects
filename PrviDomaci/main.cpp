#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <numeric>
#include <mgl2/mgl.h>
#include <mgl2/base.h>

static Eigen::Vector3d dotA = {-3, -2, 1};
static Eigen::Vector3d dotB = {1, -3, 1};
static Eigen::Vector3d dotC = {1, 3, 1};
static Eigen::Vector3d dotD = {-4, 1, 1};
static Eigen::Vector3d dotE = {4, -2, 1};
static Eigen::Vector3d dotF = {4, 1, 1};
static std::vector<Eigen::Vector3d> globalInPts = {dotA, dotB, dotE, dotF, dotC, dotD};

static Eigen::Vector3d dotAp = {-2, -2, 1};
static Eigen::Vector3d dotBp = {3, -2, 1};
static Eigen::Vector3d dotCp = {3, 2, 1};
static Eigen::Vector3d dotDp = {-2, 2, 1};
static Eigen::Vector3d dotEp = {6, -1, 1};
static Eigen::Vector3d dotFp = {6, 1, 1};
static std::vector<Eigen::Vector3d> globalOutPts = {dotAp, dotBp, dotEp, dotFp, dotCp, dotDp};

Eigen::Matrix3d naive_algorithm(Eigen::Vector3d A, Eigen::Vector3d B,
                                Eigen::Vector3d C, Eigen::Vector3d D,
                                Eigen::Vector3d Ap, Eigen::Vector3d Bp,
                                Eigen::Vector3d Cp, Eigen::Vector3d Dp){
    // matrix with column A, B, C
    Eigen::Matrix3d Mat1;
    Mat1 << A, B, C;

    // solving linear system for alpha, beta, gama
    Eigen::Vector3d abg = (Mat1.transpose() * Mat1).ldlt().solve(Mat1.transpose() * D);
    A *= abg(0);
    B *= abg(1);
    C *= abg(2);

    Eigen::Matrix3d P1;
    P1 << A, B, C;

    // matrix with column A', B', C'
    Eigen::Matrix3d Mat2;
    Mat2 << Ap, Bp, Cp;

    // solving linear system for alpha', beta', gama'
    Eigen::Vector3d apbpgp = (Mat2.transpose() * Mat2).ldlt().solve(Mat2.transpose() * Dp);
    Ap *= apbpgp(0);
    Bp *= apbpgp(1);
    Cp *= apbpgp(2);

    Eigen::Matrix3d P2;
    P2 << Ap, Bp, Cp;

 //return from function P = P2 * P1^(-1)
 return P2 * P1.inverse();
}

Eigen::MatrixXd calculate_p(Eigen::Vector3d point, Eigen::Vector3d point_p){
    Eigen::MatrixXd p(2, 9);
    p << 0, 0, 0, -point_p(2)*point(0), -point_p(2)*point(1), -point_p(2)*point(2), point_p(1)*point(0), point_p(1)*point(1), point_p(1)*point(2),
         point_p(2)*point(0), point_p(2)*point(1), point_p(2)*point(2), 0, 0, 0, -point_p(0)*point(0), -point_p(0)*point(1), -point_p(0)*point(2);
 return p;
}

Eigen::Matrix3d DLT(std::vector<Eigen::Vector3d> in_points, std::vector<Eigen::Vector3d> out_points){

    unsigned long n = in_points.size();

    //calc matrix dimn 2x9 for every correspondence Mi <-> Mi'
    std::vector<Eigen::MatrixXd> Ps(n);
    for(unsigned long i = 0; i < n; i++){
        Ps.at(i) = calculate_p(in_points.at(i), out_points.at(i));
    }

    //making matrix A
    Eigen::MatrixXd A(2*n, 9);
    for(unsigned long i = 0; i < 2*n; ){
        for(unsigned k = 0; k < 2; k++,i++){
            for(unsigned j = 0; j < 9; j++){
                A(i,j) = Ps.at(i/2)(k,j);
            }
        }
    }

    //SVD decomposition of matrix A plus calc matrix V
    Eigen::JacobiSVD<Eigen::MatrixXd> svdA(A, Eigen::ComputeFullV);
    Eigen::MatrixXd tmpMat = svdA.matrixV();

    //taking the last column of matrix V and putting in vector tmpP
    Eigen::VectorXd tmpP(9);
    for(int i = 0; i < 9; i++){
        tmpP(i) = tmpMat(i, 8);
    }

    //making matrix P 3x3 from vector tmpP
    Eigen::Matrix3d P;
    P << tmpP(0), tmpP(1), tmpP(2), tmpP(3), tmpP(4), tmpP(5), tmpP(6), tmpP(7), tmpP(8);

 return P;
}

Eigen::Matrix3d normalize(std::vector<Eigen::Vector3d> points){

    unsigned long n = points.size();

    //making P(x,y) from (x1 : x2 : x3)
    std::vector<Eigen::Vector2d> pts(n);
    for(unsigned long i = 0; i < n; i++){
        pts.at(i)(0) = points.at(i)(0)/points.at(i)(2);
        pts.at(i)(1) = points.at(i)(1)/points.at(i)(2);
    }

    //calc center of system of pts
    double x = 0, y = 0;
    for(unsigned long i = 0; i < n; i++){
        x += pts.at(i)(0);
        y += pts.at(i)(1);
    }
    std::vector<double> center = {x/n, y/n};

    //matrix of translation G
    Eigen::Matrix3d G;
    G << 1, 0, center.at(0)*-1, 0, 1, center.at(1)*-1, 0, 0, 1;

    //adding z=1, making homogeneous coordinates
    std::vector<Eigen::Vector3d> pts3d(n);
    for(unsigned long i = 0; i < n; i++){
        pts3d.at(i)(0) = pts.at(i)(0);
        pts3d.at(i)(1) = pts.at(i)(1);
        pts3d.at(i)(2) = 1.0;
    }

    //translating all pts using G
    for(unsigned long i = 0; i < n; i++){
        pts3d.at(i) = G * pts3d.at(i);
    }

    //calc avg dist of points at the moment
    double dist = 0;
    for(unsigned long i = 0; i < n; i++){
        dist += (std::sqrt(pts3d.at(i)(0)*pts3d.at(i)(0) + pts3d.at(i)(1)*pts3d.at(i)(1)));
    }
    dist /= n;

    //scaling matrix tmpS, scaling in origin
    Eigen::Matrix3d tmpS;
    tmpS << sqrt(2)/dist, 0, 0, 0, sqrt(2)/dist, 0, 0, 0, 1;

    //scaling matrix S, scaling in any point
    Eigen::Matrix3d S = G.inverse() * tmpS * G;

    //normalization matrix
    Eigen::Matrix3d T = S * G;

 return T;
}

Eigen::Matrix3d nDLT(std::vector<Eigen::Vector3d> in_points, std::vector<Eigen::Vector3d> out_points){

    unsigned long n = in_points.size();

    //calc T and T', normalization matrices for inPts and outPts
    Eigen::Matrix3d T = normalize(in_points);
    Eigen::Matrix3d Tp = normalize(out_points);

    //normalizing every in and out point
    std::vector<Eigen::Vector3d> inPts(n), outPts(n);
    for(unsigned long i = 0; i < n; i++){
        inPts.at(i) = T * in_points.at(i);
        outPts.at(i) = Tp * out_points.at(i);
    }

    //calc matrix P' and from there, calc P = T'^(-1)*P'*T
    Eigen::Matrix3d Pp = DLT(inPts, outPts);
    Eigen::Matrix3d P = Tp.inverse()*Pp*T;

 return P;
}

void image(mglGraph *gr) {
    unsigned long n = globalInPts.size();
    mglData x(n), y(n), z(n);
    for(unsigned long i = 0; i < n; i++){
        x.a[i] = globalInPts.at(i)(0);
        y.a[i] = globalInPts.at(i)(1);
        z.a[i] = globalInPts.at(i)(2);
    }

    mglData x2(6), y2(6), z2(6), c2(6);
    for(unsigned long i = 0; i < n; i++){
        x2.a[i] = globalOutPts.at(i)(0);
        y2.a[i] = globalOutPts.at(i)(1);
        z2.a[i] = globalOutPts.at(i)(2);
    }

    gr->Title("DLT & nDLT");
    gr->SetRanges(0,12, 0, 12, 0, 1);
    gr->Box();
    gr->Dots(x,y,z);
    gr->Dots(x2,y2,z2);
}

int main(){
    mglGraph gr(0, 1920, 1080);
    image(&gr);
    gr.WriteFrame("6_tacaka.png");

    Eigen::Matrix3d na, dlt, ndlt;

    na = naive_algorithm(dotA, dotB, dotC, dotD, dotAp, dotBp, dotCp, dotDp);
    dlt = DLT(globalInPts, globalOutPts);
    ndlt = nDLT(globalInPts, globalOutPts);

    std::cout << "Naive algorithm:\n" << na << std::endl << std::endl;

    std::cout << "DLT algorithm:\n" << dlt << std::endl << std::endl;

    std::cout << "NDLT algorithm:\n" << ndlt << std::endl << std::endl;

 return 0;
}
