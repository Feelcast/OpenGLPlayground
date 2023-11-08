
#include "vectors.cpp"

bool areColliding(Particle p1, Particle p2){
    vec dif = p2.pos - p1.pos;
    double distsq = sqNorm(dif);
    double colRad = (p1.r + p2.r)*(p1.r + p2.r);
    bool result = distsq < 1.05*colRad;
    return result;
}

bool areColliding(vec v1, vec v2, Particle p){
    vec on = orthoNormal(v1,v2);
    vec unitary = normalized(v2-v1);
    vec middle = (v2+v1)/2.0;
    vec dvec = p.pos-middle;
    double normalProy = dot(dvec, on);
    double tProy = dot(dvec, unitary);
    double lenghtSq = sqNorm(v2-v1);
    bool result = normalProy*normalProy < (p.r*1.01)*(p.r*1.01) && tProy*tProy < 1.01*lenghtSq/4.0;
    return result;
}

int areColliding(Box b, Particle p){
    //bVertex index from 0 to 3
    BoxVertex bVertex = b.getVertex();
    int result = 1*areColliding(bVertex.vertex[0],bVertex.vertex[1],p) +
                    2*areColliding(bVertex.vertex[1],bVertex.vertex[2],p) +
                    3*areColliding(bVertex.vertex[2],bVertex.vertex[3],p) +
                    4*areColliding(bVertex.vertex[3],bVertex.vertex[0],p);
    return result;
}

void calculateCol(Particle &p1, Particle &p2){
    //I did the momentum and energy equations on the frame where particle 2 is static, momentum transfer to the second particle is on the normal
    vec normal = normalized(p2.pos - p1.pos);
    vec relv = p1.vel - p2.vel;
    double tMass = p1.mass + p2.mass;
    vec v2f = normal*(2*p1.mass/tMass)*dot(relv, normal);
    vec v1f = v2f*(-p2.mass/p1.mass);
    p1.vel = p1.vel + v1f;
    p2.vel = p2.vel + v2f;
}

void calculateCol(Box &b, Particle &p, int numcase){
    vec normal(0,0);
    BoxVertex bVertex = b.getVertex();
    switch (numcase) {
    case 1:
        normal = orthoNormal(bVertex.vertex[0],bVertex.vertex[1]);
        break;
    case 2:
        normal = orthoNormal(bVertex.vertex[1],bVertex.vertex[2]);
        break;
    case 3:
        normal = orthoNormal(bVertex.vertex[2],bVertex.vertex[3]);
        break;
    case 4:
        normal = orthoNormal(bVertex.vertex[3],bVertex.vertex[0]);
        break;
    default:
        break;
    }

    vec relv = b.vel - p.vel;
    double tMass = b.mass + p.mass;

    vec v2f = normal*(2*b.mass/tMass)*dot(relv, normal);
    vec v1f = v2f*(-p.mass/b.mass);

    if(b.mechanic){
        b.vel = b.vel + v1f;
    }
    if (p.mechanic){
        p.vel = p.vel + v2f;
    }
    
}

void particleCols(std::vector<Particle> &particles){
    int s = particles.size();
    for (int i = 0; i<s; i++){
        Particle &pi = particles[i];
        for (int j = 0; j<i; j++){
            Particle &pj = particles[j];
            if(areColliding(pi,pj)){
                calculateCol(pi,pj);
            }
        }
    }
}

void particleCols2(std::vector<Particle> &particles){
    //implementar quadtrees
    // dividir el espacio en celdas de 20x20
    int s = particles.size();
    for (int i = 0; i<s; i++){
        Particle &pi = particles[i];
    }
}

void particleBoxCols(std::vector<Particle> &particles, std::vector<Box> &boxes){
    //implementar quadtrees
    int sp = particles.size();
    int sb = boxes.size();
    for (int i = 0; i<sp; i++){
        Particle &p = particles[i];
        for (int j = 0; j<sb; j++){
            Box &b = boxes[j];
            int numcase = areColliding(b,p);
            if(numcase!=0){
                calculateCol(b,p,numcase);
            }
        }
    }
}

double forceInvSq(double dsq, double m1, double m2, double k){
    double f = 0;
    if (dsq!=0){
        f = (k*m1*m2)/(dsq);
    }
return f;
}

vec forceApply(vec v1, vec v2, double m1, double m2, int fid){
    vec v = normalized((v2-v1));
    double dsq = sqNorm(v2-v1);
    double f = 0;
    double k = 4000;
    switch(fid){
        case 1:
        f = forceInvSq(dsq,m1,m2,k);
        break;
    }
    return v*f;
}

vec forceGrav(vec v1, vec v2, double m1, double m2){
    return forceApply(v1,v2,m1,m2,1);
}

vec acSum(std::vector<Particle> &particles, int n, vec posDelta, vec force(vec x1, vec x2, double m1, double m2)){
    int s = particles.size();
    vec ac(0,0);
    Particle &pn = particles[n]; 
    for (int i=0; i<s; i++){
        Particle pi = particles[i];
        if (i!=n){
            ac = ac + force(pn.pos + posDelta, pi.pos, pn.mass, pi.mass)/pn.mass;
        }           
    }  
    if (pn.mechanic){
        pn.ac = ac;
    }
    return ac;
}

void RK4inter(std::vector<Particle> &particles, int n, double h, bool forceSim, vec force(vec x1, vec x2, double m1, double m2)){
    Particle &p1 = particles[n];
    vec zero(0,0);
    vec r1 = p1.pos;
    vec v1 = p1.vel;
    vec kr1 = v1;
    vec kv1 = acSum(particles,n,zero, force);
    vec kr2 = v1 + kv1*h/2.0;
    vec kv2 = acSum(particles,n,kr1*h/2.0, force);
    vec kr3 = v1 + kv2*h/2.0;
    vec kv3 = acSum(particles,n,kr2*h/2.0, force);
    vec kr4 = v1 + kv3*h;
    vec kv4 = acSum(particles,n,kr3*h, force);
    if (p1.mechanic){
    p1.vel = p1.vel + (kv1+kv2*2+kv3*2+kv4)*h/6.0;
    }
    p1.pos = p1.pos + (kr1+kr2*2+kr3*2+kr4)*h/6.0;
}

void constantVelSim(std::vector<Particle> &particles, int n, double h){
    Particle &p1 = particles[n];
    if(p1.vel[0]!=0 && p1.vel[1]!=0){
    p1.pos = p1.pos + p1.vel*h;
    }
}

void boxRK4(std::vector<Box> &boxes, int n, double h){
    Box &b = boxes[n];
    if(b.vel[0]!=0 && b.vel[1]!=0){
        b.pos = b.pos + b.vel*h;
    }
    
}

void particleDynamics(std::vector<Particle> &particles, bool forceSim){
    double h = 0.001;
    int s = particles.size();
    for (int i=0; i<s; i++){
        if (forceSim){
            RK4inter(particles, i, h, forceSim, forceGrav);    
        }
        else{
            constantVelSim(particles, i, h);
        }
    }
    particleCols(particles);
}

void boxDynamics(std::vector<Box> &boxes){
    double h = 0.001;
    int s = boxes.size();
    for (int i=0; i<s; i++){
        boxRK4(boxes,i,h);         
    }
}

void updateTraces(std::vector<Particle> &particles){
    for(Particle &p : particles){
        if(p.trace.size()<1000){
            p.trace.push_back(p.pos);
        }
        else{
            p.trace.erase(p.trace.begin());
            p.trace.push_back(p.pos);
        }
    }
}

void createGas(Container c, std::vector<Box> &boxes, std::vector<Particle> &particles){
const int xstep = 32;
const int ystep = 32;
int  n = c.lenght/xstep;
int  m = c.height/ystep;
vec x(1,0);
vec y(0,1);
vec ipos(0,0);

for (int  i = 0; i<n;i++){
    for (int j = 0; j<m; j++){
        ipos  = x*(-c.lenght/2.0+10)+y*(c.height/2.0-10) + x*xstep*i - y*ystep*j; 
        Particle pt;
        pt.pos = ipos;
        int vx = rand()%100 - 50;
        int vy = rand()%100 - 50;
        pt.vel = vec(vx,vy);
        pt.mass = 1;
        pt.r = 3;
        particles.push_back(pt);
    }
}

registerContainer(c,boxes);
}

void createFluid(double speed, double density){

}
// Read data from a file and return a vector of Particle objects
std::vector<Particle> readParticlesFromFile(const std::string& filename) {
    std::vector<Particle> particles; // Vector to store Particle objects

    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Error: Failed to open the file." << std::endl;
        return particles; // Return an empty vector in case of an error
    }

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream lineStream(line);
        std::vector<double> values;
        std::string cell;
        while (std::getline(lineStream, cell, ',')) {
            try {
                double value = std::stod(cell);
                values.push_back(value);
            } catch (const std::invalid_argument& e) {
                std::cerr << "Error: Invalid data format in the file." << std::endl;
                return particles; // Return an empty vector in case of an error
            }
        }
        if (values.size() != 7) {
            std::cerr << "Error: Each row should contain 7 values." << std::endl;
            return particles; // Return an empty vector in case of an error
        }
        Particle particle;
        particle.pos = vec(values[0], values[1]);
        particle.vel = vec(values[2], values[3]);
        particle.mass = values[4];
        particle.r = values[5];
        if (values[6]!=0){
            particle.mechanic = false;
        }
        particles.push_back(particle);
    }

    return particles;
}

