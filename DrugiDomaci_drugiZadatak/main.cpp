#include <iostream>
#include <GL/glut.h>
#include <Eigen/Dense>
#include <vector>
#include <QQuaternion>

#define TIMER_ID 0
#define TIMER_INTERVAL 10
static double hours;
static bool animation_active;
static int timer_active;

static void on_timer(int value) {
    if (value != 0)
        return;

    hours += 1;
    if(hours > 100)
        return;

    glutPostRedisplay();

    if (timer_active)
        glutTimerFunc(50, on_timer, 0);
}

static void on_keyboard(unsigned char key, int, int) {
    switch (key) {
    case 'g':
    case 'G':
        if (!timer_active) {
            glutTimerFunc(50, on_timer, 0);
            timer_active = 1;
            animation_active = true;
        }
        break;

    case 's':
    case 'S':
        timer_active = 0;
        break;

    case 27:
        exit(0);
    }
}

static void on_reshape(int width, int height) {
    glViewport(0, 0, width, height);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(60, width*1.0 / height, 1, 1500);
}

static void draw_coosys(float scale) {
    glBegin(GL_LINES);
    glColor3f (1, 0, 0);
    glVertex3f(scale, 0, 0);
    glVertex3f(0, 0, 0);

    glColor3f (0, 1, 0);
    glVertex3f(0, scale, 0);
    glVertex3f(0, 0, 0);

    glColor3f (0, 0, 1);
    glVertex3f(0, 0, scale);
    glVertex3f(0, 0, 0);
    glEnd();
}

static void draw_object() {
    glColor3f(1, 1, 1);
    glutWireTetrahedron();
}

typedef std::pair<Eigen::Vector3d, double> returnType; //shorter

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

    Eigen::Vector3d p;
    p << A(2,1)-A(1,2), A(0,2)-A(2,0), A(1,0)-A(0,1);

    double teta = std::asin(p.norm()/2.0);

    p = p.normalized();

 return returnType(p.transpose(),teta);
}

Eigen::Matrix3d Rodrigez(Eigen::Vector3d p, double teta){

    Eigen::Matrix3d px;
    px <<  0.0,  -p(2),  p(1),
           p(2),  0.0,  -p(0),
          -p(1),  p(0),  0.0;

    Eigen::Matrix3d E;
    E << 1.0, 0.0, 0.0,
         0.0, 1.0, 0.0,
         0.0, 0.0, 1.0;

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
    p = p.normalized();
    Eigen::Vector3d xyz = std::sin(teta/2) * p;

 return Eigen::Vector4d{xyz(0), xyz(1), xyz(2), w};
}

returnType Q2AxisAngle(Eigen::Vector4d q){

    q = q.normalized();
    double teta = 2*std::acos(q(3)); //2*arccos(w)
    Eigen::Vector3d v;
    Eigen::Vector3d tmp{q(0), q(1), q(2)};
    if(std::abs(q(3)) == 1.0)
        v << 1.0, 0.0, 0.0;
    else
        v = tmp.normalized();

 return returnType(v.transpose(), teta);
}

static void on_display(void) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(5, 5, 5, 0, 0, 0, 0, 1, 0);

    glPushMatrix();
    glTranslatef(1, 1, 3);
    Eigen::Matrix3d A = Euler2A(3*M_PI/4, M_PI/2, 2*M_PI/3);
    returnType aa1 = axisAngle(A);
    Eigen::Vector4d q1 = AxisAngle2Q(aa1.first, aa1.second);
    glRotated(aa1.second * 180 / M_PI, aa1.first(0), aa1.first(1), aa1.first(2));
    draw_coosys(1);
    draw_object();
    glPopMatrix();

    //Pocetak
    //Pozicija (1, 1, 3)
    //Ojlerovi uglovi: 3*Pi/4, Pi/2, 2*Pi/3
    //Kraj
    //Pozicija (3, 2, 0)
    //Ojlerovi uglovi: Pi/2, Pi/3, -Pi/4


    glPushMatrix();
    glTranslatef(3, 2, 0);
    Eigen::Matrix3d A2 = Euler2A(M_PI/2, M_PI/3, -M_PI/4);
    returnType aa2 = axisAngle(A2);
    Eigen::Vector4d q2 = AxisAngle2Q(aa2.first, aa2.second);
    glRotated(aa2.second * 180 / M_PI, aa2.first(0), aa2.first(1), aa2.first(2));
    draw_coosys(1);
    draw_object();
    glPopMatrix();

    glPushMatrix();
    draw_coosys(10);
    glPopMatrix();

    glPushMatrix();
    if(animation_active) {
        Eigen::Vector3d C = (1 - hours / 100) * Eigen::Vector3d(1, 1, 3) + hours / 100 * Eigen::Vector3d(3, 2, 0);
        glTranslated(C(0), C(1), C(2));
        Eigen::Vector4d q = Slerp(q1, q2, 100, hours);
        returnType aa = Q2AxisAngle(q);
        glRotated(aa.second * 180 / M_PI, aa.first(0), aa.first(1), aa.first(2));
        draw_coosys(1);
        draw_object();
    }
    glPopMatrix();

    glutSwapBuffers();
}

int main(int argc, char *argv[]) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);

    glutInitWindowSize(800, 600);
    glutInitWindowPosition(100, 100);
    glutCreateWindow(argv[0]);

    glutKeyboardFunc(on_keyboard);
    glutReshapeFunc(on_reshape);
    glutDisplayFunc(on_display);

    glClearColor(0.25f, 0.25f, 0.25f, 0);
    glEnable(GL_DEPTH_TEST);

    hours = -1;
    timer_active = 0;
    animation_active = false;

    glutMainLoop();

 return 0;
}
