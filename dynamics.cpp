
#include "vectors.cpp"

double forceInvSq(double d, double m1, double m2, double k){
    double f = 0;
    if (d!=0){
        f = (k*m1*m2)/(d*d);
    }
return f;
}

vec forceApply(Particle& p1, Particle& p2, int fid){
    vec v1 = p1.pos;
    vec v2 = p2.pos;
    double m1 = p1.mass;
    double m2 = p2.mass;
    vec v = normalized((v2-v1));
    double d = distance(v1,v2);
    double f = 0;
    double k = 400;
    switch(fid){
        case 1:
        f = forceInvSq(d,m1,m2,k);
        break;
    }
    return v*f;
}

vec forceApply(vec v1, vec v2, double m1, double m2, int fid){
    vec v = normalized((v2-v1));
    double d = distance(v1,v2);
    double f = 0;
    double k = 400;
    switch(fid){
        case 1:
        f = forceInvSq(d,m1,m2,k);
        break;
    }
    return v*f;
}

vec forceGrav(vec v1, vec v2, double m1, double m2){
    return forceApply(v1,v2,m1,m2,1);
}

vec vnext(Particle p2,Particle p1,vec f(vec x2,vec x1, double m2, double m1),double h){
  vec x2 = p2.pos;
  vec x1 = p1.pos;
  double m2 = p2.mass;
  double m1 = p1.mass;
    vec k1= f(x2,x1,m2,m1);
    vec k2 = f(x2+ k1*h/2.,x1,m2,m1);
    vec k3 = f(x2+ k2*h/2.,x1,m2,m1);
    vec k4 = f(x2+ k2*h,x1,m2,m1);
    return x2 + (k1+k2*2+k3*2+k4)*h/6.;
}

vec xnext(Particle p,double h){
  vec v = p.vel;
    vec k1 = v;
    vec k2 = v + k1*h/2.;
    vec k3 = v + k2*h/2.;
    vec k4 = v + k3*h;
    return v + (k1+k2*2+k3*2+k4)*h/6.;
}

vec nextPosGrav(Particle &p,Particle other,double h){
  vec v = vnext(p,other,forceGrav,h);
  p.vel=v;
  return xnext(p,h);
}

void update(Particle &p2,Particle &p1,double h){
p1.pos= nextPosGrav(p1,p2,h);
p2.pos= nextPosGrav(p2,p1,h);
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
        if (values.size() != 5) {
            std::cerr << "Error: Each row should contain 5 values." << std::endl;
            return particles; // Return an empty vector in case of an error
        }
        Particle particle;
        particle.pos = vec(values[0], values[1]);
        particle.vel = vec(values[2], values[3]);
        particle.mass = values[4];
        particles.push_back(particle);
    }

    return particles;
}

