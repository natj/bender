// main interface for C++ Bender

#ifndef CPPE_H
#define CPPE_H


#include <iostream>
#include <string>
#include <vector>


#include <o2scl/exception.h>
#include <o2scl/cli.h>
#include <o2scl/hdf_file.h>

#include "ns.h"


namespace cppe {

    class cppe_class {

        protected:

        // error handler for each thread
        o2scl::err_hnd_cpp error_handler;

        // arguments sent to the command-line
        std::vector<std::string> run_args;

        // command-line interface
        o2scl::cli cl;

        // main otuput stream
        std::ofstream scr_out;

        // main ns parameters needed for ns-object
        double mass, rad, spin, incld;


        // parameter object for the 'set' command
        o2scl::cli::parameter_double p_mass;
        o2scl::cli::parameter_double p_rad;
        o2scl::cli::parameter_double p_spin;
        o2scl::cli::parameter_double p_incl;


        // set up the 'cli' object
        void setup_cli();

        // nstar nss;
        // metric ms;

        int init(std::vector<std::string> &sv, bool itive_com);


        public:

            // main class containing everything
            cppe_class();

            int rtrace(nstar &ns, metric &m);

            // command-line argument parser
            void run(int argc, char *argv[]);


    };


}



#endif

