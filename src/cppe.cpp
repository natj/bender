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
    std::cout.precision(7);


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


    return;
}


