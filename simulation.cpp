#include "dynamics.cpp"
#include <GL/freeglut.h>
// k= 500
//cgs system
// global
//Quadtree qt(-100,100,800,800,1000);
std::vector<Particle> particles;
std::vector<Box> boxes;
std::vector<std::vector<vec>> partPositions;
std::vector<std::vector<vec>> boxPositions;
std::vector<std::vector<vec>> boxVels;
std::vector<std::vector<vec>> boxForces;
std::vector<LightRay> rays;
MediumMatrix opticalObjects;
long unsigned int time_int = 0;
bool isDragging = false;
vec dragStart;
int draggedObject;
//config
bool traceFlag = true;
bool forceSim = false;
bool rts = false;
bool optics = false;
//simulation constants
int frame = 0;
int frameLimit = 6100;
double h = 0.001;
//window size
const double xsc = 1280;
const double ysc = 720;
// button
bool play = false;

/*
void equiSystem(double m1, double m2, double m3, double r,double k){
vec zero(0,0);
double M = m1 + m2 + m3;
double raux =  r*(1.5);
double L = sqrt(3.0)*r;
double omega = sqrt(k*M/(L*L*L));

vec r1 = fromPolar(r,0);
double r1cm = (m2+m3)*raux/M;
vec v1 = normalized(rotationAntiClockW(r1))*r1cm*omega;

vec r2 = fromPolar(r,120);
double r2cm = (m1+m3)*raux/M;
vec v2 = normalized(rotationAntiClockW(r2))*r2cm*omega;

vec r3 = fromPolar(r,240);
double r3cm = (m1+m2)*raux/M;
vec v3 = normalized(rotationAntiClockW(r3))*r3cm*omega;

particles.push_back({r1,v1,zero,m1});
particles.push_back({r2,v2,zero,m2});
particles.push_back({r3,v3,zero,m3});
}
*/
void updateDynamics(){
    particleDynamics(particles, forceSim,h);
    particleBoxCols(particles, boxes);
    boxDynamics(boxes,h);
}

void readPositions(int frameNumber){
    int ps = particles.size();
    int bs  = boxes.size();
    for (int i = 0; i<ps;i++){
        particles[i].pos = partPositions[frameNumber][i];
    }
    /*
    for (int i = 0; i<bs;i++){
        boxes[i].pos = boxPositions[frameNumber][i];
    }
    */
}

void circle (vec pos, double cr, double r, double g, double b){
    const int n = 100;
    const double PI = 3.141592654;
    double px  = pos[0];
    double py = pos[1];
	glBegin(GL_TRIANGLE_FAN);
    glColor3f(r,g,b);
	glVertex2f(px,py);
	for ( int i=0;i<n+1;i++){
	    glVertex2f(px+cr*cosf(i*2*PI/n),py+cr*sinf(i*2*PI/n));
	}
    glColor3f(1.0,1.0,1.0);
	glEnd();
}

void circleUnfilled (vec pos, double cr){
    const int n = 100;
    const double PI = 3.141592654;
    double px  = pos[0];
    double py = pos[1];
	glBegin(GL_LINE_STRIP);
	for ( int i=0;i<n+1;i++){
	    glVertex2f(px+cr*cosf(i*2*PI/n),py+cr*sinf(i*2*PI/n));
	}
	glEnd();
}
void renderOpticalObjects(){
    for(OpticCircle c : opticalObjects.circles){
        circleUnfilled(c.pos, c.radius);
    }
}

void renderRay(LightRay ray){
    glBegin(GL_LINE_STRIP);
    glColor3f(ray.getColor().r, ray.getColor().g, ray.getColor().b);
    for(MediumVertex m : ray.interfaces){
        glVertex2f(m.pos[0],m.pos[1]);
    }
    glColor3f(1.0,1.0,1.0);
    glEnd();
}

void renderAllRays(){
    for (LightRay r : rays){
        renderRay(r);
    }
}

void line(vec v1,vec v2)
{
	glBegin(GL_LINES);
    glColor3f(0.5,0.5,0.5);
	glVertex2f(v1[0],v1[1]);
	glVertex2f(v2[0],v2[1]);
    glColor3f(1,1,1);
	glEnd();
}

void renderBox(Box b){
    BoxVertex bVertex = b.getVertexC();
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
    glBegin(GL_LINE_STRIP);
    glColor3f(0,0.6,0.7);
    for (vec v: points){
        glVertex2f(v[0],v[1]);
    }
    glColor3f(1,1,1);
	glEnd();
}

void drawPartVec(Particle p){
    vec x1 = p.pos;
    if(p.ac[0] != 0 || p.ac[1] != 0){
        vec xa = p.pos + normalized(p.ac)*5;
        colourLine(x1,xa);
    }
    if(p.vel[0] != 0 || p.vel[1] != 0){
        vec xv = p.pos + normalized(p.vel)*20;
        colourLine(x1,xv);
    }
}

void grid (float l,int num){
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

void displayText(vec pos, char *text){
    char *p;
    glPushMatrix();
    glTranslatef(pos[0], pos[1], 0);
    glScalef(1.0/10.0,1.0/10.0,1);
    for (p = text; *p; p++)
      glutStrokeCharacter(GLUT_STROKE_ROMAN, *p);
    glPopMatrix();
}

void drawButton(){
    if (play) {
        glColor3f(0.0, 1.0, 0.0); // Set the color to green when the button is pressed
    } else {
        glColor3f(1.0, 0.0, 0.0); // Set the color to red when the button is not pressed
    }
    glBegin(GL_QUADS);
    glVertex2f(500, 320); // Define the four corners of the button
    glVertex2f(520, 320);
    glVertex2f(520, 340);
    glVertex2f(500, 340);
    glEnd();
    glColor3f(1.0, 1.0, 1.0);
}

void render(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
    glOrtho(-640, 640,-360, 360, 1, -1);
	glTranslatef(0.0,0.0,0.0);
    //grid
    grid(100,20);
    //object rendering
    for (Particle p: particles){
        if (p.mechanic){
            circle(p.pos, p.r,0.7,0.1,0);
        }
        else{
            circle(p.pos, p.r,0.9,0.9,0);
        }
        
        //drawPartVec(p);
        if(traceFlag){
        renderTrace(p.trace);
        }
    }
    for (Box b: boxes){
        renderBox(b);
    }
    if (optics){
        renderAllRays();
        renderOpticalObjects();
    }
    //time text
    double time_seg = time_int/1000.0;
    char* char_time = new char[10];
    std::sprintf(char_time,"%2.3f",time_seg);
    displayText(vec(500,-320),char_time);
    //time button
    drawButton();
    // buffer swapping
    glutSwapBuffers();
}

void mouseMotion(int x, int y){
    vec cursorPos(x-xsc/2.0, ysc/2.0 - y);
    if(isDragging){
        OpticCircle &c = opticalObjects.circles[draggedObject];
        c.pos = c.pos + cursorPos - dragStart;
        dragStart = cursorPos;
    }
}

void time(int t){
    if (play){
        
        if (rts){
            //for (int i = 0; i<16;i++){
                updateDynamics();
            //}  
        }
        if (!rts && frame < frameLimit){
            readPositions(frame);
            frame++;
        }
        if (frame>= frameLimit){
            play = false;
        }
        if(time_int%96 == 0){
            if (traceFlag){
                updateTraces(particles);
            }
            /*
            if (optics){
                clearPaths(rays);
                RaySimulation(rays, opticalObjects, xsc, ysc);
            }
            */
        }
        time_int+=16;
    }   
    render();
    //redisplay
    glutPostRedisplay();
    glutTimerFunc(16,time,0);
}

void idleFunction(){
    if (play){   
        if (rts){
            updateDynamics();
        }
        if (!rts && frame < frameLimit){
            readPositions(frame);
            frame++;
        }
        if (frame>= frameLimit){
            play = false;
        }
        if(time_int%96 == 0){
            if (traceFlag){
                updateTraces(particles);
            }
            if (optics){
                clearPaths(rays);
                RaySimulation(rays, opticalObjects, xsc, ysc);
            }
        }
        time_int+=16;
    }   
    if(time_int%16 == 0){
    //se llama en vsync
    render();
    //redisplay
    glutPostRedisplay();
    }
}

void onMouseClick(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
        // Check if the mouse click coordinates are within the button's area
        if (x >= 1140 && x <= 1160 && y >= 20 && y <= 40) {
            play = !play; // Toggle the boolean variable
        }
        int index = 0;
        if(optics){
            for(OpticCircle &c: opticalObjects.circles){
                vec cursorPos(x-xsc/2.0, ysc/2.0 - y);
                if(c.contains(cursorPos)){
                    isDragging = true;
                    dragStart = cursorPos;
                    draggedObject = index;
                    continue;
                }
                index++;
            }
        }
    }
    else if (button == GLUT_LEFT_BUTTON && state == GLUT_UP){
        isDragging = false;
    }
}

int main(int argc, char** argv){
    // File containing the data
    if (rts){
        std::string filename = "particle_data.txt";
        particles = readParticlesFromFile(filename);
        //optics
        /*
        rays.push_back(LightRay(vec(-635,99),0, 390));
        rays.push_back(LightRay(vec(-635,99.5),0, 430));
        rays.push_back(LightRay(vec(-635,100),0, 565));  
        rays.push_back(LightRay(vec(-635,100.5),0, 610));
        rays.push_back(LightRay(vec(-635,101),0, 690));
        opticalObjects.circles.push_back(OpticCircle(vec(0,0), 120,1.1));
        opticalObjects.circles.push_back(OpticCircle(vec(350,0), 100,1.1));
        */
    }
    else{
        //preGenerativeMK1(0.001,frameLimit);
        initObjects(particles, boxes);
        loadSimData(partPositions, boxPositions, boxVels, boxForces);
        readPositions(0);
    }
    //traceFlag = true;
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
    //init
    glutInitWindowSize(xsc, ysc);
    glutInitWindowPosition(100, 100);
    glutCreateWindow("Simulation");
    glutDisplayFunc(render);
    glutMouseFunc(onMouseClick);
    glutMotionFunc(mouseMotion);
    //glutIdleFunc(idleFunction);
    glutTimerFunc(16,time,0);
    glutMainLoop();
    return 0;
}