#include "ns.h"

#include <fstream>



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
    double o2 = (o20 + o21*X)*pow(Osb,2);
    return 1.0 + o2*pow(cos(theta),2);
}

double nstar_obl::dRf(double theta) {
    double o2 = (o20 + o21*X)*pow(Osb,2);
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
    os << "          Rp: " << ns.Rf(0.0) << endl;
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




geo::geo(double p_aa, double p_bb, metric &p_m) {

    a = p_aa;
    b = p_bb;
    m = &p_m;

    C = pow(a, 2) + pow(b, 2); 
    Lz = a*m->sini;


    timec = 0.0;
    Xob   = 0.0;
    theta = 0.0;
    phi   = 0.0;
    cosa  = 0.0;
    opz   = 0.0;

}

void geo::set_polar(double rad, double chi) {
    a = rad*std::sin(chi);
    b = rad*std::cos(chi);

    C = pow(a, 2) + pow(b, 2); 
    Lz = a*m->sini;

    return;
}

double geo::ptim(double x, double the) {
    double enu2 = pow((m->enu(x,the)),2);
    double w = m->w(x);
    // double nu2 = m->nu2(the);

    // cout << "tim enu0:" << (m->enu(x,the))/exp(nu2*pow(x,3)) << endl;
    // cout << "tim enu2:" << exp(nu2*pow(x,3)) << endl;
    // cout << "tim enu^2:" << enu2 << endl;
    // cout << "tim enu^-12:" << 1/enu2 << endl;
    // cout << "tim w: " << w << endl;
    // cout << "tim wLz: " << (-1.0 - Lz*w) << endl;

    return (-1.0 - Lz*w)/enu2;
}

tuple<double, bool> geo::prad(double x, double the) {
    double enu2 = pow((m->enu(x,the)),2);
    double w = m->w(x);
    double B = m->B(x, the);
    double ezeta = m->ezeta(x, the);

    double sq = pow(-1.0 - Lz*w, 2) - C*pow(enu2*x/B, 2);
    bool rbool = false;

    // cout << "rad sq:" << sq << endl;
    // cout << "rad sq1:" << pow(-1.0 - Lz*w, 2) << endl;

    // cout << "rad sq2:" <<  C*pow(enu2*x/B, 2) << endl;
    if (sq < 0.0) {
        sq *= -1.0;
        rbool = true;
    }

    // return sqrt(pow(-1.0 - Lz*w,2) - C*pow(enu2,2)*pow(x,2)/pow(B,2))/ezeta;

    return make_tuple(sqrt(sq)/ezeta, rbool);
}

tuple<double, bool> geo::pthe(double x, double the) {
    double enu2 = pow((m->enu(x,the)),2);
    double B = m->B(x, the);
    double ezeta = m->ezeta(x, the);

    double sq = C-pow(Lz/sin(the),2);
    bool tbool = false;

    // cout << "the sq" << sq << endl; 
    if (sq < 0.0 ) {
        sq *= -1.0;
        tbool = true;
    }

    return make_tuple(-(enu2*pow(x,2)/ezeta/B)*sqrt(sq), tbool);
}

double geo::pphi(double x, double the) {
    double enu2 = pow((m->enu(x,the)),2);
    double w = m->w(x);
    double B = m->B(x, the);

    return (w*(-1.0 - Lz*w) + pow(enu2*x,2)*Lz/pow(B*sin(the),2))/enu2;
}


geo_rk geo::rk_step(double x, double the) {

    geo_rk rk;
    double pt, pr, py, pz, dr;
    
    dr                = -1.0/pow(x, 2);
    pt                = ptim(x, the);
    tie(pr, rk.rturn) = prad(x, the);
    tie(py, rk.tturn) = pthe(x, the);
    pz                = pphi(x, the);

    rk.dt = dr*pt/pr;
    rk.dy = dr*py/pr;
    rk.dz = dr*pz/pr;

    return rk;
}

bool geo::bender(nstar &ns) {

    
    // propagation directions
    double rsign = 1.0;
    double psign;
    if (b >= 0.0) {psign = -1.0;}
    else {psign = 1.0;}

    // leapfrog parameters
    const double h = 0.002;
    const double tol = 5.0e-7;

    // bool hit = true;
    double rr = h/256;

    // series expanded initial values
    // assuming O(rr)^1
    double tm1 = 0.0;
    double ym1 = asin(m->sini);
    double zm1 = 0.0;

    double tn = tm1 + rr*C/2.0;
    double yn = ym1 - rr*b*psign;
    double zn = zm1 + rr*a/m->sini;

    // TODO clean this
    double hi = h;
    double level = pow(2.0, 7);

    double tni, yni, zni, rri, xoi, maxri = 0.0;
    double rsigni = rsign;
    double psigni = psign;

    double tp1_o, yp1_o, zp1_o;
    double errt, erry, errz;

    int counter = 0;
    double maxr = rr;
    // double Xob = 0.0;


    // std::cout << "tn " << tn << std::endl;
    // std::cout << "yn " << yn << std::endl;
    // std::cout << "zn " << zn << std::endl;

    // std::ofstream myfile;
    // myfile.open("path.csv");

    while (true) {
        double err = 1.0;
        tni = tn;
        yni = yn;
        zni = zn;
        rri = rr;
        // xoi = Xob;
        // maxri = maxr;
        psigni = psign;
        rsigni = rsign;

        double tp1, yp1, zp1 = 0.0;
        bool level_break = true;
        while (err > tol && level_break) {
            level *= 2;
            
            // reset state
            yn = yni;
            psign = psigni;
            rsign = rsigni;
            hi = h / level;

            geo_rk rk1, rk2;
            rk1 = rk_step(rr, yn);
            rk2 = rk_step(rr + rsign*hi, yn + psign*rk1.dy*hi);

            if (rk1.tturn) {psign *= -1.0;}
            if (rk1.rturn) {rsign *= -1.0;}

            // rk21 adaptive
            tp1_o = tn + hi*rk1.dt;
            yp1_o = yn + hi*rk1.dy*psign;
            zp1_o = zn + hi*rk1.dz;

            tp1 = tn + hi*0.5*(rk1.dt + rk2.dt); 
            yp1 = yn + hi*0.5*(rk1.dy + rk2.dy)*psign;
            zp1 = zn + hi*0.5*(rk1.dz + rk2.dz);

            errt = std::abs((tp1 - tp1_o)/tp1);
            erry = std::abs((yp1 - yp1_o)/yp1);
            errz = std::abs((zp1 - zp1_o)/zp1);
            err = std::max({10.0*errt, erry, errz});

            // std::cout << "tp1 " << tp1 << std::endl;
            // std::cout << "yp1 " << yp1 << std::endl;
            // std::cout << "zp1 " << zp1 << std::endl;


            if (level > 512.0) { level_break = false; }
        }

        rr += rsign*hi;

        level = max(0.5, level/4.0);

        tn = tp1;
        yn = yp1;
        zn = zp1;


        // check surface location
        // Change isoradial surface location to isotropic coordinates
        // XXX fixme
        Xob = (ns.X/ns.Rf(yn)) * m->R2iR(rr, yn);

        // std::cout << tn << ", " << yn << ", " << zn << " " << 

        //             ns.Rf(yn) << " " << Xob << " " << err << std::endl;
        // myfile << rr << ", " << tn << ", " << yn << ", " << zn << std::endl;

        // keep track of U-turns
        if (rr > maxr) { maxr = rr;}

        // terminate if we miss by 5%
        if (rr < maxr) {
            if (rr < Xob/0.95) {
                // hit = false;
                return false;
            }
        }

        if (rr >= Xob) {
            double serr = (Xob - rr)/Xob;
            counter += 1;

            if (counter > 50) {break;}

            if (abs(serr) > 1.0e-6) {
                tn = tni;
                yn = yni;
                zn = zni;
                rr = rri;
                maxr = rri;
                rsign = rsigni;
                psign = psigni;
                level *= 4.0;

            } else {
                break;
            }
        }
    }

    timec = tn;
    theta = yn;
    phi = tools::mod2pi(tools::pi - zn) - tools::pi;

    // emission angle


    // myfile.close();

    return true;
}

ostream &operator<<(ostream &os, geo &g) {

    os << "geodesic: " << endl;
    os << "    x: " << g.a << endl;
    os << "    y: " << g.b << endl;
    os << "   Lz: " << g.Lz << endl;
    os << "    C: " << g.C << endl;


    return os;
}



