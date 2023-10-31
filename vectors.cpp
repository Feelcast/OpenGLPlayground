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
        upSide = Box(pos + y*(h/2.0+2),vec(0,0),1000,4,l+8,0,false);
        downSide = Box(pos - y*(h/2.0+2),vec(0,0),1000,4,l+8,0,false);
        rightSide = Box(pos + x*(l/2.0+2),vec(0,0),1000,h,4,0,false);
        leftSide = Box(pos - x*(l/2.0+2),vec(0,0),1000,h,4,0,false);
    }
    private:

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
vec normalized(const vec& v) {
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

vec rotationClockW(vec v){
    vec vr(-v[1],v[0]);
    return vr;
}

vec orthoNormal(vec v1, vec v2){
    vec dif = v2-v1;
    return normalized(rotationClockW(dif));
}