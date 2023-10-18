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

struct Particle {
    vec pos;  // Position vector
    vec vel;  // Velocity vector
    double mass;       // Mass of the particle
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


