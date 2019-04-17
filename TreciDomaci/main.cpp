#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <mgl2/mgl.h>
#include <mgl2/base.h>

static std::vector<Eigen::Vector3d> x = {{760, 260, 1},
                                         {652, 498, 1},
                                         {752, 549, 1},
                                         {904, 424, 1},
                                         {690, 707, 1},
                                         {391, 769, 1},
                                         {689, 1034, 1},
                                         {1145, 782, 1}};

static std::vector<Eigen::Vector3d> y = {{741, 266, 1},
                                         {663, 472, 1},
                                         {717, 539, 1},
                                         {933, 474, 1},
                                         {572, 655, 1},
                                         {455, 661, 1},
                                         {559, 968, 1},
                                         {1113, 889, 1}};

/* tacke iz poruke, za proveru
static std::vector<Eigen::Vector3d> x = {{958, 38, 1},
                                         {1117, 111, 1},
                                         {874, 285, 1},
                                         {707, 218, 1},
                                         {292, 569, 1},
                                         {770, 969, 1},
                                         {770, 1465, 1},
                                         {317, 1057, 1}};
static std::vector<Eigen::Vector3d> y = {{933, 33, 1},
                                         {1027, 132, 1},
                                         {692, 223, 1},
                                         {595, 123, 1},
                                         {272, 360, 1},
                                         {432, 814, 1},
                                         {414, 1284, 1},
                                         {258, 818, 1}};*/

//Imamo zadate po 8 tacaka sa obe slike, po 4 temena 2 kvadra i to:
/* x1 -> gornje srednje teme, kutije za sibice
 * x2 -> donje levo teme
 * x3 -> donje srednje teme
 * x4 -> donje desno teme
 *
 * x5-x8 isto kao x1-x4, samo temena kutije
 *
 * y1-y8 isto kao x1-x8, samo za desnu sliku
 * Videce se sa slike ako nije najjasnije
 *
 * Sada treba izracunati preostala temena sa slika(vidljiva i nevidljiva)
 * x9  -> gornje levo teme, kutije za sibice
 * x10 -> gornje desno teme
 * x11 -> gornje nevidljivo teme
 * x12 -> donje nevidljivo teme
 *
 * respektivno za x13-x16, samo temena kutije
 *
 * respektivno za y9-y16, samo temena kutija sa desne slike
*/

void izracunajPreostaleTacke(){

    x.push_back({649, 210, 1}); //x9
    x.push_back({919, 142, 1}); //x10
    x.push_back({814, 93, 1}); //x11

    //x12, nevidljivo teme kutije za sibice
    Eigen::Vector3d xb = (x.at(2).cross(x.at(0))).cross(x.at(3).cross(x.at(9)));
    Eigen::Vector3d yb = (x.at(0).cross(x.at(9))).cross(x.at(8).cross(x.at(10)));
    Eigen::Vector3d app = x.at(10).cross(xb);
    Eigen::Vector3d bpp = x.at(1).cross(yb);
    Eigen::Vector3d XX = app.cross(bpp);
    XX /= XX(2);
    x.push_back(XX);

    x.push_back({373, 444, 1}); //x13
    x.push_back({1167, 467, 1}); //x14

    //x15, nevidljivo gornje teme velike kutije
    xb = (x.at(4).cross(x.at(13))).cross(x.at(6).cross(x.at(7)));
    yb = (x.at(4).cross(x.at(12))).cross(x.at(6).cross(x.at(5)));
    app = x.at(12).cross(xb);
    bpp = x.at(13).cross(yb);
    XX = app.cross(bpp);
    XX /= XX(2);
    x.push_back(XX);

    //x16, nevidljivo donje teme velike kutije
    xb = (x.at(6).cross(x.at(4))).cross(x.at(7).cross(x.at(13)));
    yb = (x.at(4).cross(x.at(13))).cross(x.at(12).cross(x.at(14)));
    app = x.at(14).cross(xb);
    bpp = x.at(5).cross(yb);
    XX = app.cross(bpp);
    XX /= XX(2);
    x.push_back(XX);


    y.push_back({677, 200, 1}); //y9
    y.push_back({965, 206, 1}); //y10
    y.push_back({903, 139, 1}); //y11

    //y12 nevidljivo teme kutije za sibice, desna slika
    xb = (y.at(2).cross(y.at(0))).cross(y.at(3).cross(y.at(9)));
    yb = (y.at(0).cross(y.at(9))).cross(y.at(8).cross(y.at(10)));
    app = y.at(10).cross(xb);
    bpp = y.at(1).cross(yb);
    Eigen::Vector3d YY = app.cross(bpp);
    YY /= YY(2);
    y.push_back(YY);

    y.push_back({463,357, 1}); //y13
    y.push_back({1148, 582, 1}); //y14
    y.push_back({998, 298, 1}); //y15

    //y16 nevidljivo teme velike kutije, desna slika
    xb = (y.at(6).cross(y.at(4))).cross(y.at(7).cross(y.at(13)));
    yb = (y.at(4).cross(y.at(13))).cross(y.at(12).cross(y.at(14)));
    app = y.at(14).cross(xb);
    bpp = y.at(5).cross(yb);
    YY = app.cross(bpp);
    YY /= YY(2);
    y.push_back(YY);
}

Eigen::Matrix3d normalize(const std::vector<Eigen::Vector3d> & points){

    size_t n = points.size();

    //making P(x,y) from (x1 : x2 : x3)
    /*std::vector<Eigen::Vector2d> pts(n);
    for(size_t i = 0; i < n; i++){
        pts.at(i)(0) = points.at(i)(0)/points.at(i)(2);
        pts.at(i)(1) = points.at(i)(1)/points.at(i)(2);
    }*/

    //calc center of system of pts
    double xdist = 0, ydist = 0;
    for(size_t i = 0; i < n; i++){
        xdist += points.at(i)(0);
        ydist += points.at(i)(1);
    }
    std::vector<double> center = {xdist/n, ydist/n};

    //matrix of translation G
    Eigen::Matrix3d G;
    G << 1, 0, center.at(0)*-1, 0, 1, center.at(1)*-1, 0, 0, 1;

    //adding z=1, making homogeneous coordinates
    std::vector<Eigen::Vector3d> pts3d(n);
    for(size_t i = 0; i < n; i++){
        pts3d.at(i)(0) = points.at(i)(0);
        pts3d.at(i)(1) = points.at(i)(1);
        pts3d.at(i)(2) = 1.0;
    }

    //translating all pts using G
    for(size_t i = 0; i < n; i++){
        pts3d.at(i) = G * pts3d.at(i);
    }

    //calc avg dist of points at the moment
    double dist = 0;
    for(size_t i = 0; i < n; i++){
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

Eigen::Matrix3d calculateF(){
    //matrica 8 jednacina, dobijenih iz korespodencije yT * F * x = 0
    Eigen::MatrixXd jed8(8, 9);
    for(size_t i = 0, n = 8; i < n; i++){
        Eigen::VectorXd tmp(9);
        tmp << y.at(i)(0) * x.at(i)(0),
               y.at(i)(0) * x.at(i)(1),
               y.at(i)(0) * x.at(i)(2),
               y.at(i)(1) * x.at(i)(0),
               y.at(i)(1) * x.at(i)(1),
               y.at(i)(1) * x.at(i)(2),
               y.at(i)(2) * x.at(i)(0),
               y.at(i)(2) * x.at(i)(1),
               y.at(i)(2) * x.at(i)(2);
        for(long j = 0; j < 9; j++)
            jed8(long(i), j) = tmp(j);
    }

    //svd dekompozicija matrice jed8, racunamo samo matricu V, zbog ustede vremena
    Eigen::JacobiSVD<Eigen::MatrixXd> svdJed8(jed8, Eigen::ComputeFullV);
    Eigen::MatrixXd fullV = svdJed8.matrixV();

    //uzimamo poslednju kolonu matrice V i smestamo u vector F, iz koga cemo u sledecem koraku dobiti fundamentalnu matricu
    Eigen::VectorXd FVector(9);
    for(int i = 0; i < 9; i++){
        FVector(i) = fullV(i, 8);
    }

    //pravimo matricu F i ubacujemo brojeve iz vektora F
    Eigen::Matrix3d FF;
    FF << FVector(0),
          FVector(1),
          FVector(2),
          FVector(3),
          FVector(4),
          FVector(5),
          FVector(6),
          FVector(7),
          FVector(8);

 return FF;
}

Eigen::Vector3d izracunajEpipol(const Eigen::MatrixXd & matrix){
    Eigen::Vector3d e;
    for(int i = 0; i < 3; i++)
        e(i) = matrix(i, 2);
 return e;
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> dekomponuj(const Eigen::Matrix3d & T0){
    Eigen::HouseholderQR<Eigen::Matrix3d> qr(T0);
    Eigen::MatrixXd Q = qr.householderQ();
    Eigen::MatrixXd R = qr.matrixQR().triangularView<Eigen::Upper>();
    //std::cout << "Provera dekompozicije matrice T0:\n" << Q*R << '\n' << std::endl;

    //matrica kalibracije K = R^(-1)
    R = R.inverse();
    //deljenje sa elementom [3,3], da bi matrica kalibracije bila po formuli
                                                                            /*K = [[dx,  s,  x],
                                                                                   [0,  dy,  y],
                                                                                   [0,   0,  1]*/
    //radimo ovo u main-u, nakon provere T0 = K*A
    //R /= R(2,2);

 return std::pair<Eigen::MatrixXd, Eigen::MatrixXd>(Q, R);
}

int main() {

    izracunajPreostaleTacke();

    //normalizovanje tacaka
    Eigen::Matrix3d T = normalize(x);
    Eigen::Matrix3d Tp = normalize(y);

    unsigned long n = x.size();

    /*
    std::cout << "X tacke:\n" << std::endl;
    for(unsigned long i = 0; i < n; i++){
        std::cout << x.at(i).transpose() << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Y tacke:\n" << std::endl;
    for(unsigned long i = 0; i < n; i++){
        std::cout << y.at(i).transpose() << std::endl;
    }
    std::cout << std::endl;
    */

    for(unsigned long i = 0; i < n; i++){
        x.at(i) = T * x.at(i);
        y.at(i) = Tp * y.at(i);
    }

    //F predstavlja fundamentalnu matricu
    Eigen::Matrix3d F = calculateF();

    std::cout << "Fundamentalna:\n" << F << '\n' << std::endl;

    //svd dekompozicija fundamentalne matrice F, rastavljam na UDV
    Eigen::JacobiSVD<Eigen::MatrixXd> svdF(F, Eigen::ComputeFullU | Eigen::ComputeEigenvectors | Eigen::ComputeFullV);
    //pravim dijagonalnu matricu od vectora D, tako da poslednja vrednost bude 0
    Eigen::Vector3d sv = svdF.singularValues();
    Eigen::Matrix3d svMatrix;
    svMatrix << sv(0), 0, 0, 0, sv(1), 0, 0, 0, 0;
    //matrica U
    Eigen::Matrix3d svdFfullU = svdF.matrixU();
    //matrica V
    Eigen::Matrix3d svdFfullV = svdF.matrixV();

    //levi epipol (epipol 1)
    Eigen::Vector3d e1 = izracunajEpipol(svdFfullV);
    std::cout << "e1: " << (e1/e1(2)).transpose() << '\n' << std::endl;

    //desni epipol (epipol 2)
    Eigen::Vector3d e2 = izracunajEpipol(svdFfullU);
    std::cout << "e2: " << (e2/e2(2)).transpose()  << '\n' << std::endl;

    //Nova fundamentalna matrica, F1, dobijena mnozenjem U * izmenjeno D * V^T
    //U daljem toku zadatka koristimo nju umesto matrice F
    Eigen::Matrix3d F1 = svdFfullU * svMatrix * svdFfullV.transpose();
    std::cout << "Fundamentalna 1:\n" << F1 << '\n' <<  std::endl;

    std::cout << "Determinanta fundamentalne matrice 1: " << F1.determinant() << '\n' << std::endl;

    //kanonske matrice kamera
    Eigen::MatrixXd T1(3, 4);
    T1 << Eigen::Matrix3d::Identity(), Eigen::Vector3d::Zero();
    std::cout << "T1:\n" << T1 << '\n' <<  std::endl;

    //parametri kamere 1

    //pozicija kamere 1 (C1)  T1 * C1 = (0, 0, 0)^T, resavam pomocu svd dekompozicije
    std::cout << "Pozicija kamere 1: " <<  T1.jacobiSvd(Eigen::ComputeFullV).matrixV().col(3).transpose() << '\n' << std::endl;

    //T1 = [T01 | -T02*C2] -> T01 = E
    Eigen::MatrixXd T01 = Eigen::Matrix3d::Identity();
    //dekompozicija matrice T01, na matricu rotacije A1 i matricu kalibracije K1
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> A1K1 = dekomponuj(T01);

    std::cout << "A1:\n" << A1K1.first << '\n' << std::endl;
    std::cout << "K1:\n" << A1K1.second << '\n' << std::endl;


    //matrica vektorskog mnozenja sa e2
    Eigen::Matrix3d E2;
    E2 << 0, -e2(2), e2(1),
          e2(2), 0, -e2(0),
          -e2(1), e2(0), 0;
    //matrica kamere 2 -> T2[M | e2] = T2[E2*F1 | e2] => M = E2*F1
    Eigen::MatrixXd T2(3, 4), M = E2*F1;
    T2 << M, e2;
    std::cout << "T2:\n" << T2 << '\n' << std::endl;

    //parametri kamere 2

    //pozicija kamere 2 (C2) T2 * C2 = (0, 0, 0)^T, resavam pomocu svd dekompozicije
    Eigen::VectorXd C2 = T2.jacobiSvd(Eigen::ComputeFullV).matrixV().col(3);
    //racunam 3D koordinate kamere (delim sa 4. parametrom)
    //C2 /= C2(3);
    std::cout << "Pozicija kamere 2: " <<  C2.transpose() << '\n' << std::endl;

    //T02 za kameru 2   T2 = [T02 | -T02*C2] => T02 = M
    Eigen::MatrixXd T02 = M;
    //Pseudo inverz matrice T02 (A^+ = (A^T * A)^(-1) * A^T)
    Eigen::MatrixXd T02pi = (T02.transpose() * T02).inverse() * T02.transpose();
    //dekompozicija matrice T02, na matricu rotacije A2 i matricu kalibracije K2
    std::pair<Eigen::MatrixXd, Eigen::MatrixXd> A2K2 = dekomponuj(T02pi);

    std::cout << "Provera K2*A2:\n" << A2K2.second * A2K2.first.transpose() << '\n' << std::endl;
    //matrica rotacije A = Q^(-1) = Q^(T)
    std::cout << "A2:\n" << A2K2.first.transpose() << '\n' << std::endl;
    std::cout << "K2:\n" << A2K2.second / A2K2.second(2,2) << '\n' << std::endl;


    //triangulacija
    Eigen::MatrixXd X(n,4);

    for(unsigned long i = 0; i < n; i++){
        Eigen::Matrix4d A;

        A <<  x.at(unsigned(i))(1)*T1.row(2) - T1.row(1),
             -x.at(unsigned(i))(0)*T1.row(2) + T1.row(0),
              y.at(unsigned(i))(1)*T2.row(2) - T2.row(1),
             -y.at(unsigned(i))(0)*T2.row(2) + T2.row(0);

        //Eigen::Vector4d B = Eigen::Vector4d::Zero();
        Eigen::Vector4d svd = A.jacobiSvd(Eigen::ComputeFullV).matrixV().col(A.cols()-1);
        //predstavljanje homogenih koordinata u 3d prostoru
        svd /= svd(3);
        //predstavljanje homogenih koordinata u 2d prostoru
        //svd /= svd(2);
        X.row(long(i)) = svd;
    }

    std::cout << "X za svih 16 tacaka:\n" << X << std::endl;

    //  plot
    /*mglGraph gr(0, 1000, 1000);
    mglData x_axis(n/2), y_axis(n/2), z_axis(n/2);
    mglData x_sibice_axis(n/2), y_sibice_axis(n/2), z_sibice_axis(n/2);
    long j = 0, k = 0;
    //u prvi niz se smestaju tacke kutije za kameru, a u drugi kutije sibica
    for(long i = 0; i < long(n); i++){
        if(i < 4 || (i > 7 && i < 12)){
            x_sibice_axis.a[k] = X(i, 0);
            //Obrtanje y koordinate tacke, zbog pravilnog prikaza na slici, pikseli krecu iz gornjeg-levog coska
            y_sibice_axis.a[k] = X(i, 1) * (-1);
            z_sibice_axis.a[k] = X(i, 2);
            k++;
        }
        else{
            x_axis.a[j] = X(i, 0);
            //Obrtanje y koordinate tacke, zbog pravilnog prikaza na slici, pikseli krecu iz gornjeg-levog coska
            y_axis.a[j] = X(i, 1) * (-1);
            z_axis.a[j] = X(i, 2);
            j++;
        }
    }

    gr.Title(":D");
    gr.SetRanges(768, 773, -472, -467, 1, 1);
    gr.Axis("xy");
    gr.Box();
    //iscrtavanja kutije sibica plavom bojom
    gr.Dots(x_sibice_axis,y_sibice_axis,z_sibice_axis, "b2");
    //iscrtavanja kutije za kameru crvenom bojom
    gr.Dots(x_axis,y_axis,z_axis, "r2");
    gr.WriteFrame("slika.png");*/


 return 0;
}
