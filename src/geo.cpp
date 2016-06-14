#include "geo.h"


using namespace std;
using namespace cppe;


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

// struct geo_rk {
// 
//     // double pa;
//     // double pr;
//     // double pt;
//     // double pp;
//     // double dr;
// 
//     double dt;
//     double dy;
//     double dz;
//     
//     bool tturn;
//     bool rturn;
// 
// }

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

void geo::bender(nstar ns) {

    
    // propagation directions
    double rsign = 1.0;
    double psign;
    if (b >= 0.0) {psign = -1.0;}
    else {psign = 1.0;}

    // leapfrog parameters
    const double h = 0.002;
    const double tol = 5.0e-7;

    bool hit = true;
    double rr = h/256;

    // series expanded initial values
    // assuming O(rr)^1
    double tm1 = 0.0;
    double ym1 = asin(m->sini);
    double zm1 = 0.0;

    double tn = tm1 + rr*C/2.0;
    double yn = ym1 - rr*b*psign;
    double zn = zm1 + rr*Lz;

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
            hi = h/level;

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

            errt = abs(tp1 - tp1_o);
            erry = abs(yp1 - yp1_o);
            errz = abs(zp1 - zp1_o);

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
        // Rgm, dtR = Rgmf(yn);
        // Xobi = X/Rgm;
        // Xob = Xobi*B/enu;
        Xob = ns.X/ns.Rf(yn) * m->R2iR(ns.X, yn);


        // keep track of U-turns
        if (rr > maxr) { maxr = rr;}

        // terminate if we miss by 5%
        if (rr < maxr) {
            if (rr < Xob/0.95) {
                hit = false;
                break;
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



    return;
}

ostream &operator<<(ostream &os, geo &g) {

    os << "geodesic: " << endl;
    os << "    x: " << g.a << endl;
    os << "    y: " << g.b << endl;
    os << "   Lz: " << g.Lz << endl;
    os << "    C: " << g.C << endl;


    return os;
}
