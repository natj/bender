#include "ns.h"
#include "geo.h"


using namespace std;
using namespace cppe;


geo::geo(double p_aa, double p_bb, metric &p_m) {

    a = p_aa;
    b = p_bb;
    m = &p_m;

    C = pow(a, 2) + pow(b, 2); 
    Lz = a*m->sini;

}

double geo::ptim(double x, double the) {
    double enu2 = pow((m->enu(x,the)),2);
    double w = m->w(x);
    double nu2 = m->nu2(the);

    cout << "enu0:" << (m->enu(x,the))/exp(nu2*pow(x,3)) << endl;
    cout << "enu2:" << exp(nu2*pow(x,3)) << endl;
    cout << "enu^2:" << enu2 << endl;
    cout << "w: " << w << endl;

    return (-1.0 + Lz*w)/enu2;
}

double geo::prad(double x, double the) {
    double enu2 = pow((m->enu(x,the)),2);
    double w = m->w(x);
    double B = m->B(x, the);
    double ezeta = m->ezeta(x, the);

    return sqrt(pow(-1.0 + Lz*w,2) - C*pow(enu2,2)*pow(x,2)/pow(B,2))/ezeta;
}

double geo::pthe(double x, double the) {
    double enu2 = pow((m->enu(x,the)),2);
    // double w = m->w(x);
    double B = m->B(x, the);
    double ezeta = m->ezeta(x, the);

    return (enu2*pow(x,2)/ezeta/B)*sqrt(C-pow(Lz/sin(the),2));
}

double geo::pphi(double x, double the) {
    double enu2 = pow((m->enu(x,the)),2);
    double w = m->w(x);
    double B = m->B(x, the);
    // double ezeta = m->ezeta(x, the);

    return (w*(-1.0 + Lz*w) + pow(enu2*x,2)*Lz/pow(B*sin(the),2))/enu2;
}


ostream &operator<<(ostream &os, geo &g) {

    os << "geodesic: " << endl;
    os << "    x: " << g.a << endl;
    os << "    y: " << g.b << endl;
    os << "   Lz: " << g.Lz << endl;
    os << "    C: " << g.C << endl;


    return os;
}
