#include "dynamics.cpp"
#include <GL/freeglut.h>
// k= 500

// global
std::vector<Particle> particles;
std::vector<Box> boxes;
bool traceFlag = false;
bool forceSim = false;

void equiSystem(double m1, double m2, double m3, double r,double k){
vec zero(0,0);
double M = m1 + m2 + m3;
double raux =  r*(1.5);
double L = sqrt(3.0)*r;
double omega = sqrt(k*M/(L*L*L));

vec r1 = fromPolar(r,0);
double r1cm = (m2+m3)*raux/M;
vec v1 = normalized(rotationClockW(r1))*r1cm*omega;

vec r2 = fromPolar(r,120);
double r2cm = (m1+m3)*raux/M;
vec v2 = normalized(rotationClockW(r2))*r2cm*omega;

vec r3 = fromPolar(r,240);
double r3cm = (m1+m2)*raux/M;
vec v3 = normalized(rotationClockW(r3))*r3cm*omega;

particles.push_back({r1,v1,zero,m1});
particles.push_back({r2,v2,zero,m2});
particles.push_back({r3,v3,zero,m3});
}

void updateDynamics(){
    particleDynamics(particles, forceSim);
    particleBoxCols(particles, boxes);
    boxDynamics(boxes);
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

void renderBox(Box b){
    BoxVertex bVertex = b.getVertex();
    glBegin(GL_QUADS);
    glColor3f(0.75,0.75,0);
    glVertex2f(bVertex.vertex[0][0],bVertex.vertex[0][1]);
    glVertex2f(bVertex.vertex[1][0],bVertex.vertex[1][1]);
    glVertex2f(bVertex.vertex[2][0],bVertex.vertex[2][1]);
    glVertex2f(bVertex.vertex[3][0],bVertex.vertex[3][1]);
    glColor3f(1,1,1);
    glEnd();
}

void colourLine(vec v1,vec v2)
{
	glBegin(GL_LINES);
    glColor3f(0.8,0,0);
	glVertex2f(v1[0],v1[1]);
	glVertex2f(v2[0],v2[1]);
    glColor3f(1,1,1);
	glEnd();
}

void renderTrace(std::vector<vec> points){
    glBegin(GL_LINES);

    for (vec v: points){
        glVertex2f(v[0],v[1]);
    }
	glEnd();
}

void drawPartVec(Particle p){
    vec x1 = p.pos;
    if(p.ac[0] != 0 && p.ac[1] != 0){
        vec xa = p.pos + normalized(p.ac)*5;
        colourLine(x1,xa);
    }
    if(p.vel[0] != 0 && p.vel[1] != 0){
        vec xv = p.pos + normalized(p.vel)*20;
        colourLine(x1,xv);
    }
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

//TO DO
void killFarObjects(){

}

void render(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
    glOrtho(-640, 640,-360, 360, 1, -1);
	glTranslatef(0.0,0.0,0.0);
    grid(200,10);
    for (Particle p: particles){
        circle(p.pos, p.r);
        //drawPartVec(p);
        //renderTrace(p.trace);
    }
    for (Box b: boxes){
        renderBox(b);
    }
    glutSwapBuffers();
}

void time(int t){
    updateDynamics();
    t++;
    if(t%15 == 0){
        render();
        glutPostRedisplay();
    }

    if(t%90 == 0 && traceFlag){
        updateTraces(particles);
        t = 0;
    }
    glutTimerFunc(1,time,t);
}

int main(int argc, char** argv){
    // File containing the data
    std::string filename = "particle_data.txt";
    particles = readParticlesFromFile(filename);
    traceFlag = true;
    //box creation
    boxes.push_back(Box(vec(0,0), vec(0,0),100,160,160,0,true));
    //forceSim = true;
    //equiSystem(10,10,10,40,1000);

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