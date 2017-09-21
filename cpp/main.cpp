#include <iostream>
#include "cppe.h"


/*  Main program to handle the command line arguments etc.
    Eventually all of the MPI + openMP initializations 
    should be also done here */
int main(int argc, char *argv[]){

	// suppress warning
	(void)argc; (void)argv;

    // add MPI inits/finalizes here

	std::cout << "hello from cppender" << std::endl;

    cppe::cppe_class m;
    m.run(argc, argv);

	return 0;
}
