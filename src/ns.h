#ifndef NS_H
#define NS_H

#include <iostream>
#include <cmath>
#include <tuple>
#include "tools.h"


namespace cppe {

    /* General ns class that computes different physical parameters
        from given observables */
    // XXX should this be a general unusable parent class?
    class nstar {

        // protected:
            // functions for derived params

        public:
            
            // basic ns parameters/observables
            double mass;
            double rad;
            double spin;

            // derived parameters
            //-------------------------------------------------- 

            // dimensionless angular velocity (in units of Keplerian vel)
            double Osb;   

            // Radius in Sch units
            double X;

            // dimensionless moment of inertia I
            double imom;

            // dimensionless angular momentum J/M^2
            double jmom;  

            // main class
            nstar(double p_mass, double p_rad, double p_spin);

            // normalized radius function and derivative against colatitude
            virtual double Rf(double theta); 
            virtual double dRf(double theta); 


    };
    

    // oblate ns
    class nstar_obl : public nstar {

        private:

            // AlGendy & Morsink 2014 EoS/R polynomial params
            const double o20 = -0.788;
            const double o21 = 1.030;

        public:

            nstar_obl(double p_mass, double p_rad, double p_spin);

            double Rf(double theta);
            double dRf(double theta);
    };

    std::ostream &operator<<(std::ostream &os, nstar &ns);

    // General parent metric class
    // uses formalism by Butterworth & Ipser
    class metric {

            public:

            // observer specific values
            //-------------------------------------------------- 
            
            // input inclination
            double incld, incl;

            // sin(incl)
            double sini;

            // Sch radius
            double Rg;

            // Dimensionless R
            double X;

            // image scaling factor from Sch units to cm^2
            double imgscale;

            // spacetime specific stuff
            //-------------------------------------------------- 

            // angular velocity of the local intertial frame
            double wspin; 

            // pressure and energy quadrupole moment
            double beta;  
            double quad;  

            metric(nstar &ns, double p_incld);

            // second order (HT) metric expansion terms
            virtual double nu2(double theta);
            virtual double B2(double theta);
            virtual double zeta2(double theta);

            // metric functions
            virtual double w(double rr);
            virtual double enu(double rr, double theta);
            virtual double B(double rr, double theta);
            virtual double ezeta(double rr, double theta);

            virtual double R2iR(double rr, double theta);

            
            // Geodesic of the metric
            // class geo {
            //     public:
            //         // Constants of motion for the geodesic
            //         double C;
            //         double Lz;

            //         // main class
            //         geo(double p_aa, double p_b);

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

            //         //  basic photon parameters
            //         double a;
            //         double b;

            // };
            // forward-decleration
            // class geo;

            




    };

    class metric_sch : public metric {
        public:
            metric_sch(nstar &ns, double p_incld);
    };

    std::ostream &operator<<(std::ostream &os, metric &m);


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

            void set_polar(double rad, double chi);


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

            bool bender(nstar &ns);



        private:
           
            // pointer to the metric of the geodesic
            metric* m;


    };
    std::ostream &operator<<(std::ostream &os, geo &ph);



}

#endif
