#include<iostream>
#include<cmath>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

class vec {
    public:
    // Constructors
    vec() : x(0.0), y(0.0) {}
    vec(double x, double y) : x(x), y(y) {}

    // Overloaded operators
    vec operator+(const vec& other) const {
        return vec(x + other.x, y + other.y);
    }

    vec operator-(const vec& other) const {
        return vec(x - other.x, y - other.y);
    }

    vec operator*(double scalar) const {
        return vec(x * scalar, y * scalar);
    }

    vec operator/(double scalar) const {
        if (scalar != 0.0) {
            return vec(x / scalar, y / scalar);
        } else {
            // Handle division by zero gracefully
            std::cerr << "Error: Division by zero." << std::endl;
            return vec(x, y);
        }
    }

    // Bracket initialization
    double& operator[](int index) {
        if (index == 0) {
            return x;
        } else if (index == 1) {
            return y;
        } else {
            throw std::out_of_range("Index out of range");
        }
    }

    const double& operator[](int index) const {
        if (index == 0) {
            return x;
        } else if (index == 1) {
            return y;
        } else {
            throw std::out_of_range("Index out of range");
        }
    }

    // Output operator
    friend std::ostream& operator<<(std::ostream& os, const vec& vec) {
        os << "(" << vec.x << ", " << vec.y << ")";
        return os;
    }

    // Getter for x coordinate
    double getX() const {
        return x;
    }

    // Getter for y coordinate
    double getY() const {
        return y;
    }

private:
    double x, y;
};

// Dot product as a standalone function
double dot(const vec& v1, const vec& v2) {
    return v1.getX() * v2.getX() + v1.getY() * v2.getY();
}

// Cross product as a standalone function
double cross(const vec& v1, const vec& v2) {
    return v1.getX() * v2.getY() - v1.getY() * v2.getX();
}

// Norm (magnitude) as a standalone function
double norm(const vec& v) {
    return sqrt(v.getX() * v.getX() + v.getY() * v.getY());
}

double sqNorm(const vec& v) {
    return v.getX() * v.getX() + v.getY() * v.getY();
}
// Normalized vector as a standalone function
vec normalized(const vec v) {
    double nor = norm(v);
    if (nor != 0.0) {
        return vec(v.getX() / nor, v.getY() / nor);
    } else {
        // Handle division by zero gracefully
        std::cerr << "Error: Normalization of a zero-length vector." << std::endl;
        return v;
    }
}

double distance(vec v1, vec v2){
    return norm(v2-v1);
}

vec fromPolar(double r, double t){
    const double PI = 3.141592654;
    double angle = t*PI/180.0;
    vec q(r*cos(angle),r*sin(angle));
    return q;
}

vec rotationAntiClockW(vec v){
    vec vr(-v[1],v[0]);
    return vr;
}

vec rotation(vec v, double t){
vec vr(v[0]*cos(t)-v[1]*sin(t),v[0]*sin(t)+v[1]*cos(t));
return vr;
}

vec orthoNormal(vec v1, vec v2){
    vec dif = v2-v1;
    return normalized(rotationAntiClockW(dif));
}

struct BoxVertex{
    vec vertex[4];
};

struct Particle {
    vec pos;  // Position vector
    vec vel;  // Velocity vector
    vec ac;
    double mass;       // Mass of the particle
    double r;
    bool mechanic = true;
    std::vector<vec> trace;
};

class Box {
    public:
    Box(): pos(0,0), vel(0,0) , mass(1), height(1),lenght(1),rotation(0) {}
    Box(vec p, vec v, double m, double h, double l, double rot, bool mec){
        const double PI = 3.141592654;
        pos = p;
        vel = v;
        mass = m;
        height = h;
        lenght = l;
        rotation = rot*PI/180.0;
        mechanic = mec;
    }
    vec pos;  // Position vector
    vec vel;  // Velocity vector
    vec ac;
    double mass;       // Mass of the particle
    double height;
    double lenght;
    double rotation;
    bool mechanic = true;

    BoxVertex getVertex(){
        vec xp(cos(rotation), sin(rotation));
        vec yp(-sin(rotation), cos(rotation));
        BoxVertex result;
        result.vertex[0] = pos - xp*lenght/2.0 + yp*height/2.0;
        result.vertex[1] = pos + xp*lenght/2.0 + yp*height/2.0;
        result.vertex[2] = pos + xp*lenght/2.0 - yp*height/2.0;
        result.vertex[3] = pos - xp*lenght/2.0 - yp*height/2.0;
        return result;
    } 
    private:
    
};

class OpticBox {
    public:
    OpticBox(vec p, double h, double l, double rot, double n){
        const double PI = 3.141592654;
        pos = p;
        height = h;
        lenght = l;
        rotation = rot*PI/180.0;
        refractionIndex = n;
    }
    vec pos;  // Position vector
    double height;
    double lenght;
    double rotation;
    double refractionIndex;

    BoxVertex getVertex(){
        vec xp(cos(rotation), sin(rotation));
        vec yp(-sin(rotation), cos(rotation));
        BoxVertex result;
        result.vertex[0] = pos - xp*lenght/2.0 + yp*height/2.0;
        result.vertex[1] = pos + xp*lenght/2.0 + yp*height/2.0;
        result.vertex[2] = pos + xp*lenght/2.0 - yp*height/2.0;
        result.vertex[3] = pos - xp*lenght/2.0 - yp*height/2.0;
        return result;
    } 
    private:
    
};

class OpticCircle{
    public:
    vec pos;
    double radius;
    double refractionIndex = 1.0;
    bool rayInside = false;
    OpticCircle(vec p, double r, double n): pos(p), radius(r), refractionIndex(n) {}
    bool contains(vec p){
        return sqNorm(p-pos)<= radius*radius;
    }
    private:
};

class Container{

    public:
    vec pos;
    double height;
    double lenght;
    //wall thickness = 4
    Box leftSide;
    Box rightSide;
    Box upSide;
    Box downSide;

    Container(vec p, double h, double l){
        pos = p;
        height = h;
        lenght = l;
        vec x(1,0);
        vec y(0,1);
        upSide = Box(pos + y*(h/2.0+2),vec(0,0),10000,4,l+8,0,false);
        downSide = Box(pos - y*(h/2.0+2),vec(0,0),10000,4,l+8,0,false);
        rightSide = Box(pos + x*(l/2.0+2),vec(0,0),10000,h,4,0,false);
        leftSide = Box(pos - x*(l/2.0+2),vec(0,0),10000,h,4,0,false);
    }
    private:

};

struct Cell {
    double x, y; // Cell's top-left corner coordinates
    double width, height; // Cell's dimensions
    std::vector<Particle*> particles; // Particles in the cell
};

struct MediumMatrix{
    std::vector<OpticBox> boxes;
    std::vector<OpticCircle> circles;
};

struct MediumVertex{
    vec pos;
    double angle;
};
struct rgbSet{
    double r;
    double g;
    double b;
};

class LightRay{
public:
vec pos;
double v = 1;
vec vel;
//wavelenght from 380 to 700
double wavelength;
double initialAngle;
std::vector<MediumVertex> interfaces;

LightRay(vec position, double ia, double w): pos(position), initialAngle(ia), wavelength(w) {
    vel = fromPolar(v,initialAngle);
    interfaces.push_back({pos, initialAngle});
    constructColor();
}

void checkMedium(MediumMatrix &m){
    for (OpticCircle &c : m.circles){
        if(c.contains(pos)){
            if (!inMedium && !c.rayInside){
                vec normal = normalized(pos-c.pos);
                vec veln = normalized(vel);
                double angle1 = dot(normal, veln);
                double angle2 = asin((currentN/c.refractionIndex)*sin(angle1));
                interfaces.push_back({pos, angle1-angle2});
                vel = rotation(vel, angle1-angle2);
                currentN =  c.refractionIndex;
                inMedium = true;
                c.rayInside = true;
            }
        }
        else{
            if (inMedium && c.rayInside){
                vec normal = normalized(pos-c.pos);
                vec veln = normalized(vel);
                double angle1 = dot(normal, veln);
                double angle2 = asin((currentN/1.0)*sin(angle1));
                interfaces.push_back({pos, angle1-angle2});
                vel = rotation(vel, angle1-angle2);
                currentN = 1.0;
                inMedium = false;
                c.rayInside = false;
            }
        }
    }
    //change the velocity and angle of incidence according to the refraction index
}

void move(MediumMatrix &m, double xsc, double ysc){
    checkMedium(m);
    pos = pos + vel;
    double limitx = xsc-10;
    double limity = ysc-10;
    if (pos[0]*pos[0]>limitx*limitx/4.0||pos[1]*pos[1]>limity*limity/4.0){
        interfaces.push_back({pos,0});
    }
}

rgbSet getColor(){return rayColor;}

private:
    double currentN = 1.00;
    bool inMedium = false;
    rgbSet rayColor = {1.0, 1.0, 1.0};

    void constructColor(){
    double gamma = 0.8;
    double intensityMax = 1.0;
    double factor = 0.0;
    double R = 0.0, G = 0.0, B = 0.0;

    if ((wavelength >= 380) && (wavelength < 440)) {
        R = -(wavelength - 440) / (440 - 380);
        G = 0.0;
        B = 1.0;
    } else if ((wavelength >= 440) && (wavelength < 490)) {
        R = 0.0;
        G = (wavelength - 440) / (490 - 440);
        B = 1.0;
    } else if ((wavelength >= 490) && (wavelength < 510)) {
        R = 0.0;
        G = 1.0;
        B = -(wavelength - 510) / (510 - 490);
    } else if ((wavelength >= 510) && (wavelength < 580)) {
        R = (wavelength - 510) / (580 - 510);
        G = 1.0;
        B = 0.0;
    } else if ((wavelength >= 580) && (wavelength < 645)) {
        R = 1.0;
        G = -(wavelength - 645) / (645 - 580);
        B = 0.0;
    } else if ((wavelength >= 645) && (wavelength <= 700)) {
        R = 1.0;
        G = 0.0;
        B = 0.0;
    }

    // Adjust intensity
    if ((wavelength >= 380) && (wavelength < 645)) {
        factor = 0.3 + 0.7 * (wavelength - 380) / (645 - 380);
    } else {
        factor = 1.0;
    }

    // Apply gamma correction
    R = std::pow(R * factor, gamma);
    G = std::pow(G * factor, gamma);
    B = std::pow(B * factor, gamma);

    rayColor = {R*intensityMax, G*intensityMax, B*intensityMax};

    }
};

void registerContainer(Container c, std::vector<Box> &boxes){
    boxes.push_back(c.leftSide);
    boxes.push_back(c.rightSide);
    boxes.push_back(c.upSide);
    boxes.push_back(c.downSide);
}


