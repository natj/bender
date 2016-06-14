#ifndef GEO_H
#define GEO_H


#include <iostream>
#include <cmath>
#include <tuple>

#include "ns.h"
// #include "tools.h"


namespace cppe {

    // Photon object 
    // class metric::geo {

    //     // protected:
    //         // functions for derived params

    //     public:
    //         
    //         // Constants of motion for the geodesic
    //         double C;
    //         double Lz;

    //         // main class
    //         geo(double p_aa, double p_bb, metric &p_m);

    //         // moments
    //         //-------------------------------------------------- 
    //         double ptim(double rr, double theta); 
    //         // double prad(double rr, double theta); 
    //         // double pthe(double rr, double theta); 
    //         // double pphi(double rr, double theta); 

    //         // derived Carter's constant values using 
    //         // radial and theta moments
    //         // double Crad(double rr, double theta); 
    //         // double Cthe(double rr, double theta); 

    //     private:
    //         
    //         //  basic photon parameters
    //         double a;
    //         double b;
    //         metric &m;

    // };
    // std::ostream &operator<<(std::ostream &os, metric::geo &ph);

    // struct geo_rk;
    struct geo_rk {

        // double pa;
        // double pr;
        // double pt;
        // double pp;
        // double dr;

        double dt;
        double dy;
        double dz;

        bool tturn;
        bool rturn;
    };


    class geo {
        public:
            //  basic photon parameters
            double a;
            double b;

            // Constants of motion for the geodesic
            double C;
            double Lz;

            double timec;
            double Xob;
            double theta;
            double phi;

            double cosa;
            double opz;

            // main class
            geo(double p_aa, double p_bb, metric &p_m);

            // moments
            //-------------------------------------------------- 
            // double ptim(metric &m, double rr, double theta); 
            double ptim(double x, double the); 
            std::tuple<double, bool> prad(double x, double the); 
            std::tuple<double, bool> pthe(double x, double the); 
            double pphi(double x, double the); 

            // derived Carter's constant values using 
            // radial and theta moments
            // double Crad(double rr, double theta); 
            // double Cthe(double rr, double theta); 

            geo_rk rk_step(double x, double the);

            void bender(nstar ns);



        private:
           
            // pointer to the metric of the geodesic
            metric* m;


    };
    std::ostream &operator<<(std::ostream &os, geo &ph);



}

#endif
