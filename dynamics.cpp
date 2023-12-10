
#include "vectors.cpp"
#include <random>

SquaredSpeedHistogram particleSpeedStatistics(std::vector<Particle> particles, double maxVelSquared, double binsNumber){
    std::vector<double> bins;
    std::vector<int> speedDistribution;
    for (int i = 1; i<=binsNumber;i++){
        bins.push_back(i*i*(maxVelSquared)/(binsNumber*binsNumber));
        speedDistribution.push_back(0);
    }
    for(Particle p: particles){
        double speedSquared = dot(p.vel,p.vel);
        double lastBinUpperLimit = 0;
        for(int j = 0; j<binsNumber;j++){
            if(speedSquared > lastBinUpperLimit && speedSquared <= bins[j]){
                speedDistribution[j]++;
            }
            lastBinUpperLimit = bins[j];
        }
    }
    SquaredSpeedHistogram result = {bins, speedDistribution};
    return result;
}

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

bool areColliding(vec m, vec n, vec u, double lsq, Particle p){
    vec dvec = p.pos-m;
    double normalProy = dot(dvec, n);
    double tProy = dot(dvec, u);
    bool result = normalProy*normalProy < (p.r*1.01)*(p.r*1.01) && tProy*tProy < 1.01*lsq/4.0;
    return result;
}

int areColliding(Box b, Particle p){
    //bVertex index from 0 to 3
    BoxMecLines bl = b.getLines();
    int result = 1*areColliding(bl.middlePoint[0], bl.normal[0], bl.unit[0], bl.lsq[0],p) +
                    2*areColliding(bl.middlePoint[1], bl.normal[1], bl.unit[1], bl.lsq[1],p) +
                    3*areColliding(bl.middlePoint[2], bl.normal[2], bl.unit[2], bl.lsq[2],p) +
                    4*areColliding(bl.middlePoint[3], bl.normal[3], bl.unit[3], bl.lsq[3],p);
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
    BoxMecLines bl = b.getLines();
    normal = bl.normal[numcase - 1];

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
        b.deltap = b.deltap + v1f*p.mass;
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

void particleBoxCols(std::vector<Particle> &particles, std::vector<Box> &boxes){
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

double forceGravPerturbative(double dsq, double m1, double m2, double k, double c){
    double f = 0;
    double d = sqrt(dsq);
    if (dsq!=0){
        f = (k*m1*m2)/(dsq) - c/(d*d*d);
    }
    return f;
}


vec forceApply(vec v1, vec v2, double m1, double m2, int fid){
    vec v = normalized((v2-v1));
    double dsq = sqNorm(v2-v1);
    double f = 0;
    double k = 4000;
    double c = 10000000;
    switch(fid){
        case 1:
            f = forceInvSq(dsq,m1,m2,k);
            break;
        case 2: 
            f = forceGravPerturbative(dsq,m1,m2,k,c);
            break;
    }
    return v*f;
}

vec forceGrav(vec v1, vec v2, double m1, double m2){
    return forceApply(v1,v2,m1,m2,1);
}
vec forceGravPrec(vec v1, vec v2, double m1, double m2){
    return forceApply(v1,v2,m1,m2,2);
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
    b.update();
    if(b.vel[0]!=0 || b.vel[1]!=0){
        b.pos = b.pos + b.vel*h;  
    }
}

void updateBoxForces(std::vector<Box> &boxes, double t){
    for (Box &b : boxes){
        b.calcForce(t);
    }
}

void particleDynamics(std::vector<Particle> &particles, bool forceSim, double h){
    int s = particles.size();
    for (int i=0; i<s; i++){
        if (forceSim){
            RK4inter(particles, i, h, forceSim, forceGravPrec);    
        }
        else{
            constantVelSim(particles, i, h);
            Particle &pi = particles[i];
            for (int j = 0; j<i; j++){
                Particle &pj = particles[j];
                if(areColliding(pi,pj)){
                    calculateCol(pi,pj);
                }
            }
        }
    }
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

double uniformDistributionSpeed(int meanSpeed){
    double randomUnit = (double) rand()/RAND_MAX;
    return randomUnit*meanSpeed*2;
}
double uniformDistributionAngle(){
    double PI = 3.141592654;
    return (rand()%100)*2*PI/100.0;
}

double maxwellBoltzmannSpeed(double meanSpeed, double mass){
     // Constants
    double energy = 0.5*mass*meanSpeed*meanSpeed;
    const double sqrt2Pi = std::sqrt(2.0 * M_PI);

    // Calculate the standard deviation
    double sigma = std::sqrt(energy / mass);

    // Create a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(0.0, sigma);

    // Generate a random speed following Maxwell-Boltzmann distribution
    double speed;
    do {
        speed = std::abs(distribution(gen));
    } while (distribution(gen) < std::exp(-mass * speed * speed / (2.0 * energy)));

    return speed;
}

void createGas(double l, double h, std::vector<Particle> &particles, double meanSpeed){
    const int xstep = 28;
    const int ystep = 28;
    int  n = l/xstep;
    int  m = h/ystep;
    vec x(1,0);
    vec y(0,1);
    vec ipos(0,0);
    for (int  i = 0; i<n;i++){
        for (int j = 0; j<m; j++){
            ipos  = x*(-l/2.0+10)+y*(h/2.0-10) + x*xstep*i - y*ystep*j; 
            Particle pt;
            pt.pos = ipos;
            double angle = uniformDistributionAngle();
            int vx = maxwellBoltzmannSpeed(meanSpeed,1.0)*cos(angle);
            int vy = maxwellBoltzmannSpeed(meanSpeed,1.0)*sin(angle);
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
            box.update();
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

void saveHistograms(std::vector<SquaredSpeedHistogram>& data, const std::string& countFilename, const std::string& binsFilename){
    std::ofstream countOutFile(countFilename);
    if (!countOutFile.is_open()) {
        std::cerr << "Error opening file: " << countFilename << std::endl;
        return;
    }
    for (const SquaredSpeedHistogram& h : data) {
        for (int count : h.speedDistribution){
            countOutFile << count <<  ",";
        }
        countOutFile << std::endl;
    }
    countOutFile.close();

    std::ofstream binsOutFile(binsFilename);
    if (!binsOutFile.is_open()) {
        std::cerr << "Error opening file: " << binsFilename << std::endl;
        return;
    }
    for (const SquaredSpeedHistogram& h : data) {
        for (double binLimit : h.bins){
            binsOutFile << binLimit <<  ",";
        }
        binsOutFile << std::endl;
    }
    binsOutFile.close();
}

void updateDynamics(std::vector<Particle>& particles, std::vector<Box>& boxes, double h, bool forceSim){
    particleDynamics(particles, forceSim, h);
    particleBoxCols(particles, boxes);
    boxDynamics(boxes,h);
}


void preGenerativeMK1(double h, int frames){
std::vector<SquaredSpeedHistogram> speedHistograms;
std::vector<Particle> particles;
std::vector<Box> boxes;
std::vector<std::vector<vec>> partPositions;
std::vector<std::vector<vec>> boxPositions;
std::vector<std::vector<vec>> boxVels;
std::vector<std::vector<vec>> boxForces;
// here the setup
bool forceSim = false;
createGas(1280,580, particles, 50);
Box upper(vec(0,290),vec(0,0),1000,4,1280,0,false);
Box lower(vec(0,-290),vec(0,0),1000,4,1280,0,false);
//Box right_wall(vec(600,0),vec(0,0),1000,576,4,0,false);
//Box left_wall(vec(-600,0),vec(0,0),1000,576,4,0,false);
Box mobile(vec(-790,0),vec(200,0),1000, 200, 200,0,true);
boxes.push_back(upper);
boxes.push_back(lower);
//boxes.push_back(right_wall);
//boxes.push_back(left_wall);
boxes.push_back(mobile);
saveParticleProperties(particles, "particle_properties.txt");
saveBoxProperties(boxes,"box_properties.txt");

for (int i = 0; i < frames*16; i++){
    if(i%16 == 0){
    std::cout<<"Generating frame: "<<i/16<<std::endl;
    std::vector<vec> tempPart;
    std::vector<vec> tempBox;
    std::vector<vec> tempBoxVels;
    std::vector<vec> tempBoxForce;
    for (const Particle& particle : particles) {
        tempPart.push_back(particle.pos);
    }
    for (const Box& box : boxes) {
        tempBox.push_back(box.pos);
        tempBoxVels.push_back(box.vel);
        tempBoxForce.push_back(box.force);
    }
    partPositions.push_back(tempPart);
    boxPositions.push_back(tempBox);
    boxVels.push_back(tempBoxVels);
    boxForces.push_back(tempBoxForce);  
    }
    if(i%160==0){
        updateBoxForces(boxes, 0.016);
        speedHistograms.push_back(particleSpeedStatistics(particles, 20000,40));
    }
    updateDynamics(particles,boxes,h, forceSim);
}
saveVectorsToFile(partPositions,"sim_particles.txt");
saveVectorsToFile(boxPositions,"sim_boxes.txt");
saveVectorsToFile(boxForces,"sim_boxes_forces.txt");
saveVectorsToFile(boxVels,"sim_boxes_vel.txt");
saveHistograms(speedHistograms,"speed_distributions.txt", "distribution_bins.txt");
}

void initObjects(std::vector<Particle>& particles, std::vector<Box>& boxes){
loadParticleProperties(particles,"particle_properties.txt");
loadBoxProperties(boxes, "box_properties.txt");
}

void loadSimData(std::vector<std::vector<vec>>& partPos,std::vector<std::vector<vec>>& boxPos,std::vector<std::vector<vec>>& boxVels,std::vector<std::vector<vec>>& boxForces){
loadFileToVectors(partPos,"sim_particles.txt");
loadFileToVectors(boxPos,"sim_boxes.txt");
//loadFileToVectors(boxVels,"sim_boxes_vel.txt");
//loadFileToVectors(boxForces,"sim_boxes_forces.txt");
}