// Basic header library containing tools, constants etc. for C++ codes
#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <numeric>
#include <vector>
#include <algorithm>
#include <cmath>

namespace tools{

    // mathematical constants
    //-------------------------------------------------- 
    const double pi=acos(-1.0);
    const double twopi= 2.0*acos(-1.0);

    // Riemann Zeta function values
    const double zeta2 = 1.6449340668482273;
    const double zeta3 = 1.2020569031595951;
    const double zeta4 = 1.0823232337111395;


    // physical constants (cgs)
    //-------------------------------------------------- 
    const double me      = 9.10938291e-28;
    const double mp      = 1.672621777e-24;
    const double k       = 1.3806488e-16;
    const double G       = 6.67384e-8;
    const double mu      = 1.32712440018e26;
    const double h       = 6.62606957e-27;
    const double eV      = 1.602176565e-12;
    const double sigmaSB = 5.67039934436e-5;
    const double sigmaTh = 6.65245853542e-25;
    const double km      = 1.0e5;

    // astrophysical constants (cgs)
    //-------------------------------------------------- 
    const double c    = 2.99792458e10;
    const double AU   = 1.49597870700e13;
    const double pc   = 3.08567758135e18;
    // const double Msun = 1.9884e33; // Measured
    const double Msun = 1.9885469e33; // Derived from mu = G*Msun

    const double constBB = 2.0*h*std::pow(1.0e3*eV/h,4)/std::pow(c,2);

    // common unit changes
    //-------------------------------------------------- 
    const double erg2kev = 1.0/(1.0e3*eV);
    const double km2kpc  = km/(1.0e3*pc);
    const double kpc2cm  = 1.0e3*pc;



    inline double deg2rad(double ang) {
        return ang*pi/180.0;
    };

    /* wrap x -> [0,max) */
    inline double wrapMax(double x, double max)
    {
        /* integer math: `(max + x % max) % max` */
        return std::fmod(max + fmod(x, max), max);
    };

    /* wrap x -> [min,max) */
    inline double wrapMinMax(double x, double min, double max)
    {
        return min + wrapMax(x - min, max - min);
    };

    inline double modpi(double x) {
        return wrapMinMax(x, -pi, pi);
    };

    // modulo 2pi using while
    // this is fastest when x has few rotations only
    inline double mod2pi(double x) {
        while (x > twopi) {x -= twopi;}
        while (x < 0.0) {x += twopi;}
        return x;
    };

    //-------------------------------------------------- 

    // Write CSV files
    /*
    std::ofstream opfile;
    opfile.open("spec.csv");

    vals[0] = 0.5;
    vals[1] = 14.1;
    
    
    std::cout << "0.0 ,";
    opfile<< "0.0 ,";
    for (j = 0; j < Ngradg2; j++)
    {
        std::cout << gradg_grid[j];
        opfile << gradg_grid[j]; 
        if (j != Ngradg2-1)
        {
            std::cout << ", ";
            opfile << ", ";
        }
    }
    std::cout << std::endl;
    opfile << std::endl;
     
    for (i = 0; i < Nenerg2; i++)
    {
       vals[3] = energ_grid[i];
        for (j=0; j < Ngradg2; j++)
        {

            if (j==0)
            {
                std::cout << energ_grid[i] << ", ";
                opfile << energ_grid[i] << ", ";
            }

            vals[2] = gradg_grid[j];
            tmp = flx.interp_linear(vals);
            // tmp = flx.interpolate(vals);

            std::cout << tmp;
            opfile << tmp;
         
            if (j != Ngradg2-1)
            {
                std::cout << ", ";
                opfile << ", ";
            }
        }
        std::cout << std::endl;
        opfile << std::endl;
    }
    opfile.close();
    */

    // Permutation sort
    template <typename T>
        std::vector<std::size_t> sort_permutation(
                const std::vector<T>& vec)
        {
            std::vector<std::size_t> p(vec.size());
            std::iota(p.begin(), p.end(), 0);
            std::sort(p.begin(), p.end(),
                    [&](std::size_t i, std::size_t j){ return vec[i] < vec[j]; });
            return p;
        }

    template <typename T>
        std::vector<T> apply_permutation(
                const std::vector<T>& vec,
                const std::vector<std::size_t>& p)
        {
            std::vector<T> sorted_vec(p.size());
            std::transform(p.begin(), p.end(), sorted_vec.begin(),
                    [&](std::size_t i){ return vec[i]; });
            return sorted_vec;
        }




 }


 #endif
