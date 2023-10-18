#include "dynamics.cpp"
#include <GL/freeglut.h>


std::vector<Particle> particles;

void updateDynamics(){

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
vec vi1(i*l,ysc);
vec vf1(i*l,0);;
line(vi1,vf1);
}
	for(int i=0;i<num;i++)
{
vec vi2(xsc,i*l);
vec vf2(0,i*l);
line(vi2,vf2);
}

}

void render(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
    glOrtho(0, 1280,0, 720, 1, -1);
	glTranslatef(0.0,0.0,0.0);

    vec p(400,400);
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

}

int main(int argc, char** argv){
    // File containing the data
    std::string filename = "particle_data.txt";
    particles = readParticlesFromFile(filename);

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