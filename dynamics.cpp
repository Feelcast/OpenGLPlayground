
#include "vectors.cpp"

bool areColliding(Particle p1, Particle p2){
    vec dif = p2.pos - p1.pos;
    double distsq = sqNorm(dif);
    double colRad = (p1.r + p2.r)*(p1.r + p2.r);
    bool result = distsq < 1.05*colRad;
    return result;
}

bool areColliding(Particle* p1, Particle* p2){
    vec dif = p2->pos - p1->pos;
    double distsq = sqNorm(dif);
    double colRad = (p1->r + p2->r)*(p1->r + p2->r);
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

class Quadtree {
public:
    Quadtree(double x, double y, double width, double height, int maxParticlesPerCell);
    ~Quadtree();
    void clear();
    void insert(Particle* particle);
    void remove(Particle* particle);
    void checkCollisions(std::vector<Particle> &particles);

private:
    void split();
    bool contains(const Cell& cell, const Particle* particle);
    void checkCollisionsInCell(const Cell& cell, std::vector<Particle> &particles);

    double x;
    double y;
    double width;
    double height;
    int maxParticlesPerCell;
    std::vector<Particle*> particles;
    Cell* cells[4];
};

Quadtree::Quadtree(double x, double y, double width, double height, int maxParticlesPerCell)
    : x(x), y(y), width(width), height(height), maxParticlesPerCell(maxParticlesPerCell) {
    for (int i = 0; i < 4; ++i) {
        cells[i] = nullptr;
    }
}

Quadtree::~Quadtree() {
    clear();
}

void Quadtree::clear() {
    particles.clear();
    for (int i = 0; i < 4; ++i) {
        if (cells[i] != nullptr) {
            delete cells[i];
            cells[i] = nullptr;
        }
    }
}

void Quadtree::insert(Particle* particle) {
    if (particles.size() < maxParticlesPerCell) {
        particles.push_back(particle);
        return;
    }

    if (cells[0] == nullptr) {
        split();
    }

    for (int i = 0; i < 4; ++i) {
        if (contains(*cells[i], particle)) {
            insert(particle);
        }
    }
}

void Quadtree::split() {
    double subWidth = width / 2;
    double subHeight = height / 2;
    cells[0] = new Cell{ x, y, subWidth, subHeight };
    cells[1] = new Cell{ x + subWidth, y, subWidth, subHeight };
    cells[2] = new Cell{ x, y + subHeight, subWidth, subHeight };
    cells[3] = new Cell{ x + subWidth, y + subHeight, subWidth, subHeight };

    // Reinsert particles into sub-cells
    for (auto& particle : particles) {
        for (int i = 0; i < 4; ++i) {
            if (contains(*cells[i], particle)) {
                insert(particle);
            }
        }
    }

    particles.clear();
}

bool Quadtree::contains(const Cell& cell, const Particle* particle) {
    return (particle->pos[0] >= cell.x && particle->pos[0] < cell.x + cell.width &&
            particle->pos[1] >= cell.y && particle->pos[1] < cell.y + cell.height);
}

void Quadtree::checkCollisions(std::vector<Particle> &particles) {
    for (auto& particle : particles) {
        for (int i = 0; i < 4; ++i) {
            if (cells[i] != nullptr) {
                checkCollisionsInCell(*cells[i], particles);
            }
        }
    }
}

void Quadtree::checkCollisionsInCell(const Cell& cell, std::vector<Particle> &particles) {
    for (size_t i = 0; i < cell.particles.size(); ++i) {
        for (size_t j = i + 1; j < cell.particles.size(); ++j) {
            Particle* p1 = cell.particles[i];
            Particle* p2 = cell.particles[j];
            if(areColliding(p1,p2)){
                //calculateColPtr(p1,p2, particles);
            }
            // Check for collision between p1 and p2
            // Implement your collision detection logic here
            // If a collision is detected, handle it as needed
            // (e.g., call calculateCol(p1, p2))
        }
    }
}

void Quadtree::remove(Particle* particle) {
    // Find the cell that contains the particle and remove it from that cell
    for (size_t i = 0; i < particles.size(); ++i) {
        if (particles[i] == particle) {
            particles.erase(particles.begin() + i);
            return;
        }
    }

    // If the particle was not found in the current cell, check sub-cells
    for (int i = 0; i < 4; ++i) {
        if (cells[i] != nullptr && contains(*cells[i], particle)) {
            remove(particle);
        }
    }
}

void RaySimulation(std::vector<LightRay> &rays, MediumMatrix &m, double xsc, double ysc){
    for (LightRay &l: rays){
        l.initialize();
        while (l.pos[0]*l.pos[0]<xsc*xsc/4.0 && l.pos[1]*l.pos[1]<ysc*ysc/4.0){
            l.move(m, xsc ,ysc);
        }
        
    }
}
void clearPaths(std::vector<LightRay> &rays){
    for (LightRay &l: rays){
        l.interfaces.clear();
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
    if(p1.vel[0]!=0 || p1.vel[1]!=0){
        p1.pos = p1.pos + p1.vel*h;
    }
}

void boxRK4(std::vector<Box> &boxes, int n, double h){
    Box &b = boxes[n];
    if(b.vel[0]!=0 || b.vel[1]!=0){
        b.pos = b.pos + b.vel*h;
    }
}

void particleDynamics(std::vector<Particle> &particles, bool forceSim, double h){
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

void boxDynamics(std::vector<Box> &boxes, double h){
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
}

void createGas(double l, double h, std::vector<Particle> &particles){
const int xstep = 28;
const int ystep = 28;
int  n = l/xstep;
int  m = h/ystep;
vec x(1,0);
vec y(0,1);
vec ipos(0,0);

    for (int  i = 0; i<n;i++){
        for (int j = 0; j<m; j++){
            ipos  = x*(-l/2.0+14)+y*(h/2.0-14) + x*xstep*i - y*ystep*j; 
            Particle pt;
            pt.pos = ipos;
            int vx = rand()%40 - 20;
            int vy = rand()%40 - 20;
            pt.vel = vec(vx,vy);
            pt.mass = 1;
            pt.r = 3;
            particles.push_back(pt);
        }
    }
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

void saveParticleProperties(const std::vector<Particle>& particles, const std::string& filename) {
    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Write particle properties to the file, one particle per line
    for (const Particle& particle : particles) {
        outFile << particle.mass << " " << particle.r << std::endl;
    }

    outFile.close();
}

void loadParticleProperties(std::vector<Particle>& particles, const std::string& filename) {
    std::ifstream inFile(filename);

    if (!inFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    particles.clear(); // Clear existing particles

    std::string line;
    while (std::getline(inFile, line)) {
        std::istringstream iss(line);
        Particle particle;

        // Parse mass and radius from the line
        if (iss >> particle.mass >> particle.r) {
            particles.push_back(particle);
        } else {
            std::cerr << "Error parsing line: " << line << std::endl;
        }
    }

    inFile.close();
}

void saveBoxProperties(const std::vector<Box>& boxes, const std::string& filename) {
    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Write box properties to the file, one box per line
    for (const Box& box : boxes) {
        outFile << box.mass << " " << box.height << " " << box.lenght << std::endl;
    }

    outFile.close();
}

void loadBoxProperties(std::vector<Box>& boxes, const std::string& filename) {
    std::ifstream inFile(filename);

    if (!inFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    boxes.clear(); // Clear existing boxes

    std::string line;
    while (std::getline(inFile, line)) {
        std::istringstream iss(line);
        Box box;

        // Parse mass, height, and length from the line
        if (iss >> box.mass >> box.height >> box.lenght) {
            boxes.push_back(box);
        } else {
            std::cerr << "Error parsing line: " << line << std::endl;
        }
    }

    inFile.close();
}

void saveVectorsToFile(const std::vector<std::vector<vec>>& data, const std::string& filename) {
    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Iterate over each vector of arrays
    for (const auto& innerVector : data) {
        // Iterate over each array in the inner vector
        for (const auto& arrayElement : innerVector) {
            // Output array elements separated by a space
            outFile << arrayElement[0] << " "<< arrayElement[1] << "|";;
        }

        outFile << "\n"; // Start a new line after each inner vector
    }

    outFile.close();
}
void loadFileToVectors(std::vector<std::vector<vec>>& data, const std::string& filename) {
    std::ifstream inFile(filename);

    if (!inFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }

    // Clear existing data
    data.clear();

    std::string line;
    while (std::getline(inFile, line)) {
        std::vector<vec> innerVector;
        std::istringstream iss(line);

        vec v;

        // Iterate over each Vec object in the line
        while (iss >> v[0] >> v[1]) {
            innerVector.push_back(v);

            // Check for the separator character '|'
            char separator;
            if (!(iss >> separator && separator == '|')) {
                std::cerr << "Error: Unexpected format in line: " << line << std::endl;
                return;
            }
        }

        data.push_back(innerVector);
    }
    inFile.close();
}

void updateDynamics(std::vector<Particle>& particles, std::vector<Box>& boxes, double h, bool forceSim){
    particleDynamics(particles, forceSim, h);
    particleBoxCols(particles, boxes);
    boxDynamics(boxes,h);
}


void preGenerativeMK1(double h, int frames){

std::vector<Particle> particles;
std::vector<Box> boxes;
std::vector<std::vector<vec>> partPositions;
std::vector<std::vector<vec>> boxPositions;

// here the setup
bool forceSim = false;
//Container m(vec(0,0),500,1000);
createGas(1280,580, particles);
Box upper(vec(0,290),vec(0,0),1000,4,1280,0,false);
Box lower(vec(0,-290),vec(0,0),1000,4,1280,0,false);
Box mobile(vec(-790,0),vec(300,0),1000, 200, 200,0,true);
boxes.push_back(upper);
boxes.push_back(lower);
boxes.push_back(mobile);

saveParticleProperties(particles, "particle_properties.txt");
saveBoxProperties(boxes,"box_properties.txt");

for (int i = 0; i < frames*15; i++){
    if(i%15 == 0){
    std::cout<<"Generating frame: "<<i/15<<std::endl;
    std::vector<vec> tempPart;
    std::vector<vec> tempBox;
    for (const Particle& particle : particles) {
        tempPart.push_back(particle.pos);
    }

    for (const Box& box : boxes) {
        tempBox.push_back(box.pos);
    }
    partPositions.push_back(tempPart);
    boxPositions.push_back(tempBox);
    }
    updateDynamics(particles,boxes,h, forceSim);
}

saveVectorsToFile(partPositions,"sim_particles.txt");
saveVectorsToFile(boxPositions,"sim_boxes.txt");
}

void initObjects(std::vector<Particle>& particles, std::vector<Box>& boxes){
loadParticleProperties(particles,"particle_properties.txt");
loadBoxProperties(boxes, "box_properties.txt");
}

void loadSimData(std::vector<std::vector<vec>>& partPos,std::vector<std::vector<vec>>& boxPos){
loadFileToVectors(partPos,"sim_particles.txt");
loadFileToVectors(boxPos,"sim_boxes.txt");
}