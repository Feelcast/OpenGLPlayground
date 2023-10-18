#include "dynamics.cpp"
#include <GL/freeglut.h>
// k= 500

std::vector<Particle> particles;

void equiSystem(double m1, double m2, double m3, double r,double k){
double M = m1 + m2 + m3;
double raux =  r*(1.5);
double L = sqrt(3.0)*r;
double omega = sqrt(k*M/(L*L*L))/1.3;

vec r1 = fromPolar(r,0);
double r1cm = (m2+m3)*raux/M;
vec v1 = normalized(rotationClockW(r1))*r1cm*omega;

vec r2 = fromPolar(r,120);
double r2cm = (m1+m3)*raux/M;
vec v2 = normalized(rotationClockW(r2))*r2cm*omega;

vec r3 = fromPolar(r,240);
double r3cm = (m1+m2)*raux/M;
vec v3 = normalized(rotationClockW(r3))*r3cm*omega;

particles.push_back({r1,v1,m1});
particles.push_back({r2,v2,m2});
particles.push_back({r3,v3,m3});
}


void updateDynamics(){
double h = 0.01;
int s = particles.size();
for (int i=0; i<s; i++){
    for (int j=0; j<i; j++){
        Particle &pi = particles[i];
        Particle &pj = particles[j];
        update(pi,pj,h);
    }
}

}

void circle (vec pos, double cr){
    const int n = 100;
    const double PI = 3.141592654;
    double px  = pos[0];
    double py = pos[1];
	glBegin(GL_TRIANGLE_FAN);
	glVertex2f(px,py);
	for ( int i=0;i<n+1;i++){
	glVertex2f(px+cr*cosf(i*2*PI/n),py+cr*sinf(i*2*PI/n));
	}
	glEnd();
}

void line(vec v1,vec v2)
{

	glBegin(GL_LINES);
	glVertex2f(v1[0],v1[1]);
	glVertex2f(v2[0],v2[1]);
	glEnd();
}

void grid (float l,int num){
double xsc = 1280;
double ysc = 720;
	for(int i=0;i<num;i++)
{
vec vi1(i*l,ysc/2.0);
vec vi1m(-i*l,ysc/2.0);
vec vf1(i*l,- ysc/2.0);
vec vf1m(-i*l,- ysc/2.0);
line(vi1,vf1);
line(vi1m,vf1m);
}
	for(int i=0;i<num;i++)
{
vec vi2(xsc/2.0,i*l);
vec vf2(-xsc/2.0,i*l);
vec vi2m(xsc/2.0,-i*l);
vec vf2m(-xsc/2.0,-i*l);
line(vi2,vf2);
line(vi2m,vf2m);
}

}

void render(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
    glOrtho(-640, 640,-360, 360, 1, -1);
	glTranslatef(0.0,0.0,0.0);
    for (Particle p: particles){
        circle(p.pos, p.mass);
    }
    grid(50,40);
    glutSwapBuffers();
}

void time(int v){
long unsigned int t = 0;
updateDynamics();
t++;
if(t%25 == 0){
render();
t = 0;
}
glutPostRedisplay();
glutTimerFunc(1,time,0);
}

int main(int argc, char** argv){
    // File containing the data
    std::string filename = "particle_data.txt";
    particles = readParticlesFromFile(filename);
    equiSystem(10,14,18,100,1000);

    //GL
    
    glutInit(&argc, argv);

    glutSetOption(GLUT_MULTISAMPLE, 8);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_MULTISAMPLE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH, GL_NICEST);
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_POINT_SMOOTH, GL_NICEST);

    glutInitWindowSize(1280, 720);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Simulation");
    glutDisplayFunc(render);
    glutTimerFunc(1,time,0);
    glutMainLoop();
    return 0;
}