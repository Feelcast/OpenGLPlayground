#include "vectors.cpp"


double forceInvSq(double d, double m1, double m2, double k){
    double f = 0;
    if (d!=0){
        f = (k*m1*m2)/(d*d);
    }
return f;
}

vec forceApply(vec v1, vec v2, double m1, double m2, int fid){
    vec v = normalized((v2-v1));
    double d = distance(v2,v1);
    double f = 0;
    double k = 10;
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

vec acSum(std::vector<Particle> &particles, int n, vec posDelta){
    int s = particles.size();
    vec ac(0,0);
    Particle &pn = particles[n]; 
    for (int i=0; i<s; i++){
        Particle pi = particles[i];
        if (i!=n){
            ac = ac + forceGrav(pn.pos + posDelta, pi.pos, pn.mass, pi.mass);
        }           
    }  
    return ac;
}

void RK4inter(std::vector<Particle> &particles, int n, double h){
    Particle &p1 = particles[n];
    vec zero(0,0);
    vec r1 = p1.pos;
    vec v1 = p1.vel;
    vec kr1 = v1;
    vec kv1 = acSum(particles,n,zero);
    vec kr2 = v1 + kv1*h/2.0;
    vec kv2 = acSum(particles,n,kr1*h/2.0);
    vec kr3 = v1 + kv2*h/2.0;
    vec kv3 = acSum(particles,n,kr2*h/2.0);
    vec kr4 = v1 + kv3*h;
    vec kv4 = acSum(particles,n,kr3*h);
    p1.vel = p1.vel + (kv1+kv2*2+kv3*2+kv4)*h/6.0;
    p1.pos = p1.pos + (kr1+kr2*2+kr3*2+kr4)*h/6.0;
}

