#include "tools.h"
#include "ns.h"

using namespace std;
using namespace cppe;

// spherical star
//-------------------------------------------------- 
nstar::nstar(double p_mass, double p_rad, double p_spin) {

    mass  = p_mass;
    rad   = p_rad;
    spin  = p_spin;

    /* NOTE: mu = G*Msun i.e. heliocentric gravitational constant
       known to much better accuracy than G or Msun alone. */
    X = tools::mu*mass/rad/tools::km/pow(tools::c,2);
    Osb = 2*tools::pi*spin*sqrt(pow(tools::km*rad,3)/(tools::mu*mass));
    
    // imom = 0.0;
    // jmom = 0.0;  

    // (AlGendy & Morsink 2014)
    imom = sqrt(X)*(1.136 - 2.53*X + 5.6*pow(X,2));
    jmom = imom*Osb*sqrt(pow(tools::c,2)*rad*tools::km/(tools::mu*mass)); 
}

double nstar::Rf(double theta) {
    return 1.0;
}
double nstar::dRf(double theta) {
    return 0.0;
}




// oblate star
//-------------------------------------------------- 

nstar_obl::nstar_obl(double p_mass, double p_rad, double p_spin) 
                : nstar(p_mass, p_rad, p_spin) { 

}

double nstar_obl::Rf(double theta) {
    static const double o2 = (o20 + o21*X)*pow(Osb,2);
    return 1.0 + o2*pow(cos(theta),2);
}

double nstar_obl::dRf(double theta) {
    static const double o2 = (o20 + o21*X)*pow(Osb,2);
    return -2.0*o2*sin(theta)*cos(theta);
}


ostream &cppe::operator<<(ostream &os, nstar &ns) {

    os << "neutron star: " << endl;
    os << " mass (Msun): " << ns.mass << endl;
    os << "    rad (km): " << ns.rad << endl;
    os << "   spin (Hz): " << ns.spin << endl;
    // os << "  incl (deg): " << ns.incld << endl;
    os << " -------------" << endl;
    os << "           X: " << ns.X << endl;
    os << "         Osb: " << ns.Osb << endl;
    os << " -------------" << endl;
    os << "        jmom: " << ns.jmom << endl;
    // os << "       wspin: " << ns.wspin << endl;
    // os << "        beta: " << ns.beta << endl;
    // os << "        quad: " << ns.quad << endl;

    return os;
}


metric::metric(nstar &ns, double p_incld) {

    incld = p_incld;
    incl = tools::deg2rad(incld);
    sini = sin(incl);

    Rg = 1.0;
    X = ns.X;

    imgscale = pow(tools::mu*ns.mass/pow(tools::c,2),2);

    wspin = 2.0*ns.jmom;

    // (AlGendy & Morsink 2014) 
    beta = 0.4454*pow(ns.Osb,2)*ns.X;
    quad = -0.11*pow(ns.Osb/ns.X,2);

}

double metric::w(double rr) {
    return wspin*pow(rr,3)*(1.0 - 3.0*rr);
}

double metric::nu2(double theta) {
    return beta/3.0 - quad*0.5*(3*pow(cos(theta),2)-1.0);
}

double metric::B2(double theta) {
    return beta;
}

double metric::zeta2(double theta) {
    return beta*(3.0*0.5*(3.0*pow(cos(theta),2)-1)/4.0 - 1.0/3.0);
}

double metric::enu(double rr, double theta) {
    return (1.0-rr/2.0)/(1.0+rr/2.0)*exp(nu2(theta)*pow(rr,3));
}

double metric::B(double rr, double theta) {
    return (1.0-rr/2.0)*(1.0+rr/2.0) + B2(theta)*pow(rr,2);
}

double metric::ezeta(double rr, double theta) {
    return (1.0-rr/2.0)*(1.0+rr/2.0)*exp(zeta2(theta)*pow(rr,2));
}

double metric::R2iR(double rr, double theta) {
    return B(rr, theta)/enu(rr, theta);
}


metric_sch::metric_sch(nstar &ns, double p_incld) 
        : metric(ns, p_incld) {

    // force spacetime into first order in rotation
    wspin = 0.0;
    beta  = 0.0;
    quad  = 0.0;
}


ostream &cppe::operator<<(ostream &os, metric &m) {

    os << "metric:       " << endl;
    os << "  incl (deg): " << m.incld << endl;
    os << "          Rg: " << m.Rg << endl;
    os << " -------------" << endl;
    os << "       wspin: " << m.wspin << endl;
    os << "        quad: " << m.quad << endl;
    os << "        beta: " << m.beta << endl;

    return os;
}
