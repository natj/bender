#include "cppe.h"


using namespace std;
using namespace cppe;


cppe_class::cppe_class() {

    // basic default ns parameters
    mass = 1.6;
    rad = 12.0;
    spin = 400.0;
    incld = 60.0;

}



void cppe_class::setup_cli() {

    // set options
    static const int nopt=1;

    o2scl::comm_option_s options[nopt]={
        {'r',"init","Set up files for writing",
            1,1,"<runname>",((string)"Initialize files."),
            new o2scl::comm_option_mfptr<cppe_class>(this,
                    &cppe_class::init),
            o2scl::cli::comm_option_both}

    };

    cl.set_comm_option_vec(nopt,options);


    p_mass.d = &mass;
    p_mass.help = "Mass of the NS (in Msun)";
    cl.par_list.insert(std::make_pair("mass", &p_mass));

    p_rad.d = &rad;
    p_rad.help = "Radius of the NS (in km)";
    cl.par_list.insert(std::make_pair("rad", &p_rad));

    p_spin.d = &spin;
    p_spin.help = "Spin of the NS (in Hz)";
    cl.par_list.insert(std::make_pair("spin", &p_spin));

    p_incl.d = &incld;
    p_incl.help = "Inclination of the NS observer (in deg)";
    cl.par_list.insert(std::make_pair("incl", &p_incl));

    std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);
    // std::cout.precision(7);
    std::cout.precision(2);


    return;
}



int cppe_class::init(std::vector<std::string> &sv, bool itive_com){

    // user-specified filename prefix
    if (sv.size()<2) {
        std::cout << "No filename given in cppe::cppe()." << std::endl;
        return o2scl::exc_efailed;
    }
    string fname_prefix=sv[1];

    // update filename with processor rank
    // MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &mpi_nprocs);
    // mpi_start_time=MPI_Wtime();
    // mpi_start_time=time(0);
    // fname_prefix+=((string)"_")+o2scl::itos(mpi_rank);

    // open main output file
    scr_out.open((fname_prefix+"_scr").c_str());
    scr_out.setf(std::ios::scientific);

    std::cout << "init()" << std::endl;

    scr_out.close();

    return 0;
}

int cppe_class::rtrace(nstar &ns, metric &m) {

    std::cout << "rtrace()" << std::endl;
    std::cout << ns;
    std::cout << m;
    std::cout << m.enu(0.2, 0.4) << std::endl;
    std::cout << m.B(0.2, 0.4) << std::endl;
    std::cout << m.ezeta(0.2, 0.4) << std::endl;
    std::cout << m.R2iR(0.2, 0.4) << std::endl;

    std::cout << "---------" << std::endl;
    geo pho(2.5, 3.0, m);
    std::cout << "C" << pho.C << std::endl;
    std::cout << "Lz" << pho.Lz << std::endl;
    // std::cout << "tim " << pho.ptim(0.2, 0.4) << std::endl;
    // std::cout << "rad " << pho.prad(0.2, 0.4) << std::endl;
    // std::cout << "the " << pho.pthe(0.2, 0.4) << std::endl;
    // std::cout << "phi " << pho.pphi(0.2, 0.4) << std::endl;


    bool hit;
    hit = pho.bender(ns);
    std::cout << "tim  " << pho.timec << std::endl;
    std::cout << "rad  " << pho.Xob << std::endl;
    std::cout << "the  " << pho.theta << std::endl;
    std::cout << "phi  " << pho.phi << std::endl;

    std::cout << "cosa " << pho.cosa << std::endl;
    std::cout << "1+z  " << pho.opz << std::endl;

    // --------------------------------------------------
    // Number of radial and angle points
    const int Nrad = 20;
    const int Nchi = 20;

    // limits for the initial image; scaled later on to fit the star
    double rmin = 0.0;
    double rmax = 12.0;
    
    const double dchi_edge = 0.01;
    const double chimin = 0.0 - dchi_edge;
    const double chimax = tools::twopi + dchi_edge;

    static const int Nedge = 31;
    

    // get radius limits for the image
    // std::vector<double> chis(Nedge), rlims(Nedge); 
    // o2scl::vector_grid(o2scl::uniform_grid_end<double>(0.0, tools::twopi, Nedge-1), chis);

    ubvector chis(Nedge), rlims(Nedge); 
    o2scl::vector_grid(o2scl::uniform_grid_end<double>(0.0, tools::twopi, Nedge-1), chis);

    static const int Nbi = 10;
    for(size_t i=0; i<Nedge; i++) {
        double rmini = rmin;
        double rmaxi = rmax;
        double rmid = 0.0;
        
        int N = 1;
        while(N<=Nbi) {
            rmid = (rmini + rmaxi)*0.5;
            
            bool hit;
            double xx = rmid*std::sin(chis[i]);
            double yy = rmid*std::cos(chis[i]);

            geo pho(xx, yy, m);
            hit = pho.bender(ns);
            
            if (hit) {
                rmini = rmid;
            } else {
                rmaxi = rmid;
            }
            N += 1;
        }
        rlims[i] = rmid;
    }

    for(size_t i=0; i<Nedge; i++){
        std::cout << i << " " << chis[i]/tools::twopi << " " << rlims[i] << std::endl;
    }

    rmax = *std::max_element(rlims.begin(), rlims.end())*1.005;
    // rmax = 8.0; // XXX debug addition
    
    // create edge function to get the exact shape of the outline
    o2scl::interp_vec<> edge_interp(Nedge, chis, rlims, o2scl::itp_linear); 

    std::cout << "interp:" << std::endl;
    for(size_t i=0; i<Nedge; i++){
        std::cout << i << " " << edge_interp.eval(chis[i]) << " " << rlims[i] << std::endl;
    }


    ubvector chi_grid(Nchi), rad_grid(Nrad);
    o2scl::vector_grid(o2scl::uniform_grid_end<double>(rmin, rmax, Nrad-1), rad_grid);  
    o2scl::vector_grid(o2scl::uniform_grid_end<double>(chimin, chimax, Nchi-1), chi_grid);  

    // o2scl::tensor_grid<ubvector, ubvector_size_t> img;
    size_t dims[3] = {Nrad, Nchi, 8};
    img.resize(3, dims);

    std::vector<double> grid;
    for(size_t i=0; i<Nrad; i++) { grid.push_back(rad_grid[i]); }
    for(size_t i=0; i<Nchi; i++) { grid.push_back(chi_grid[i]); }

    // different images
    grid.push_back(0.0);
    grid.push_back(1.0);
    grid.push_back(2.0);
    grid.push_back(3.0);
    grid.push_back(4.0);
    grid.push_back(5.0);
    grid.push_back(6.0);
    grid.push_back(7.0);

    img.set_grid_packed(grid);


    // fill with zero
    for(size_t i=0; i<Nchi; i++) {
        for (size_t j=0; j<Nrad; j++) {
            dims[0] = j;
            dims[1] = i;
            for (size_t k=0; k<8; k++) {
                dims[2] = k;
                img.set(dims, 0.0);
            }
        }
    }
            

    std::cout << "computing image.." << std::endl;
    boost::timer btimer;

    // geo pho(0.0, 0.0, m);
    for(size_t i=0; i<Nchi; i++) {
        double chi = chi_grid[i];

        if(i%10 == 0) {
            std::cout << chi << std::endl;
        }

        for (size_t j=0; j<Nrad; j++) {
            double rad = rad_grid[j];

            bool hit;

            pho.set_polar(rad, chi);

            hit = pho.bender(ns);
            
            dims[0] = j;
            dims[1] = i;

            dims[2] = 0;
            img.set(dims, pho.timec);
            dims[2] = 1;
            img.set(dims, std::sin(pho.phi));
            dims[2] = 2;
            img.set(dims, std::cos(pho.phi));
            dims[2] = 3;
            img.set(dims, pho.theta);
            dims[2] = 4;
            img.set(dims, pho.Xob);
            dims[2] = 5;
            img.set(dims, 1.0);
            dims[2] = 6;
            img.set(dims, pho.cosa);
            dims[2] = 7;
            img.set(dims, pho.opz);


            if(!hit) {
                std::cout << "break at r:" << rad << std::endl;
                break;
            }
        }
    }
    double etime = btimer.elapsed();
    std::cout << "elapsed time: " << etime << std::endl;


    




    return 0;
}


// int cppe_class::pulse(int argc, char *argv[]) {
int cppe_class::pulse(nstar &ns) {

    static const int Nt = 32;

    // time grid
    double stop_time = 1.0/ns.spin;
    ubvector times(Nt), phase(Nt);
    o2scl::vector_grid(o2scl::uniform_grid_end<double>(0.0, stop_time, Nt-1), times);
    for (std::size_t k=0; k<Nt; k++) { phase[k] = times[k] * ns.spin; };

    
    // dense cartesian image
    double xmin = -10.0;
    double xmax = 10.0;
    double ymin = -10.0;
    double ymax = 10.0;

    int Nx_dense = 30;
    int Ny_dense = 30;
    ubvector x_grid_d(Nx_dense), y_grid_d(Nx_dense);
    o2scl::vector_grid(o2scl::uniform_grid_end<double>(xmin, xmax, Nx_dense-1), x_grid_d); 
    o2scl::vector_grid(o2scl::uniform_grid_end<double>(ymin, ymax, Ny_dense-1), y_grid_d); 





    // locate star edges in cartesian image
    std::cout << "locating cartesian limits..." << std::endl;
    std::vector<int> x1s(Ny_dense), x2s(Ny_dense);
    int y1s = 0;
    int y2s = 0;
    bool img_top = false;


    // interpolation specific stuff
    img.set_interp_type(o2scl::itp_linear);

    for(size_t i=0; i<20; i++) {
        for (size_t j=0; j<20; j++) {
            size_t dims[3] = {j, i, 5};
            // dims[0] = j;
            // dims[1] = i;
            // dims[2] = 4;
            double tmp = img.get(dims);
            std::cout << "(" << i << "," << j << "; " << tmp <<") ";
        }
        std::cout << std::endl;
    }


    for (std::size_t j=0; j<Ny_dense; j++) {
        double y = y_grid_d[j];
        for (std::size_t i=0; i<Nx_dense; i++) {
            double x = x_grid_d[i];

            double rad = std::sqrt(pow(x,2) + pow(y,2));
            // double chi = tools::mod2pi(std::atan2(x,y));
            double chi = tools::mod2pi(tools::pi/2.0 - std::atan2(y,x));

            double vals[3];
            vals[0] = rad;
            // vals[1] = chi;
            vals[1] = 0.0;
            vals[2] = 5;
            double tmp = img.interp_linear(vals);

            std::cout << "rad:" << rad << " chi:" << chi << " f:" << tmp
            << " x:" << x << " y:" << y << std::endl;
            // std::cout << "(" << i << "," << j << "; " << tmp <<") ";
        }
        std::cout << std::endl;
    }



    for (std::size_t j=0; j<Ny_dense; j++) {
        double y = y_grid_d[j];

        bool left = false;
        for (std::size_t i=0; i<Nx_dense; i++) {
            double x = x_grid_d[i];

            double rad = std::sqrt(pow(x,2) + pow(y,2));
            double chi = tools::mod2pi(tools::pi/2.0 - std::atan2(y,x));

            double vals[3];
            vals[0] = rad;
            vals[1] = chi;
            vals[2] = 5;
            double hitd = img.interp_linear(vals);

            if (hitd > 0.55) {
                if (!left) { 
                    x1s[j] = i;
                    left = true;
                } else { x2s[j] = i;}

                if (!img_top) {
                    y1s = j;
                    img_top = true;
                }
                else {y2s = j;}
            }
        }
    }


    std::cout << "top:" << y1s << " " << y_grid_d[y1s] << std::endl;
    std::cout << "bot:" << y2s << " " << y_grid_d[y2s] << std::endl;
    for(std::size_t i=y1s; i<y2s; i++) {
        std::cout << "  " << i << " left:" << x1s[i] << " right:" << x2s[i] << std::endl;
    }


    for(std::size_t k=0; k<Nt; k++) {

        double t = times[k];
        std::cout << "k:" << k << " t:" << t << std::endl;

        double frame_y2, frame_y1, frame_x1, frame_x2;
        frame_y2 = y_grid_d[0];
        frame_y1 = y_grid_d[Nx_dense-1];
        frame_x1 = x_grid_d[Ny_dense-1];
        frame_y2 = x_grid_d[0];

        bool located_spot = false;

        for(std::size_t j=y1s; j<y2s; j++) {
            double y = y_grid_d[j];

            for(std::size_t i=x1s[j]; i<x2s[j]; i++) {
                double x = x_grid_d[i];

                double rad = std::sqrt(pow(x,2) + pow(y,2));
                double chi = tools::mod2pi(tools::pi/2.0 - std::atan2(y,x));

                // double edge_rad = edge_interp.eval(chi)
                double edge_rad = 6.5; // XXX debug addition

                if (rad <= edge_rad) {

                    // trace back to star
                    // phi, theta, time dtau
                    std::cout << "inside star" << std::endl;





                } // end of if inside edge
            } // end of x loop (i)
        } // end of y loop (j)
    } // end of time loop (k)


    return 0;
}



void cppe_class::run(int argc, char *argv[]) {

    o2scl::err_hnd=&error_handler;

    // process command-line arguments and run
    setup_cli();

    for (int i=0; i<argc; i++){
        run_args.push_back(argv[i]);
    }

    cl.prompt="cppe> ";
    cl.run_auto(argc, argv);

    // static spherical case
    // nstar nss(mass, rad, spin);
    // metric_sch ms(nss, incld);

    // std::cout << nss;
    // std::cout << ms;

    // std::cout << " ----- " << std::endl;
    // nss(1.4, 10.0, 600.0);
    // ms(nss, incld);

    // std::cout << nss;
    // std::cout << ms;
    // std::cout << " ----- " << std::endl;

    // rotating oblate case
    // nstar_obl ns(mass, rad, spin);
    // metric m(ns, incld);

    // std::cout << ns;
    // std::cout << m;
    
    // std::cout << m.enu(0.2, 0.4) << std::endl;
    // std::cout << m.B(0.2, 0.4) << std::endl;
    // std::cout << m.ezeta(0.2, 0.4) << std::endl;
    // std::cout << m.R2iR(0.2, 0.4) << std::endl;
    // rtrace(ns, m);
    
    std::cout << " ----- " << std::endl;
    nstar_obl nso(mass, rad, spin);
    metric mo(nso, incld);
    std::cout << nso << std::endl;
    rtrace(nso, mo);
    pulse(nso);






    return;
}


