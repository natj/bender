COMP=osx_gnu


# define the C compiler to use
ifeq ($(COMP),osx_clang)
CXX = clang++
CFLAGS = -Wall -g -stdlib=libc++ -std=c++0x
endif

ifeq ($(COMP),osx_gnu)
CXX = g++-5
CFLAGS = -Wall -Wextra -g -std=c++11 -O2
endif

ifeq ($(COMP),osx_mpi)
CXX = mpic++
CFLAGS = -Wall -Wextra -g -std=c++11 -O2 -DMCMC_MPI
endif

ifeq ($(COMP),linux_gnu)
CXX = g++
CFLAGS = -Wall -Wextra -g -std=gnu++11 -O2
endif

ifeq ($(COMP),linux_mpi)
CXX = mpic++
CFLAGS = -Wall -Wextra -g -std=gnu++11 -O2
endif

ifeq ($(COMP),cluster_gnu)
CXX = mpic++
CFLAGS = -Wall -Wextra -g -std=gnu++11 -O2 -I/home/jatnat/lib/include -L/home/jatnat/lib/lib
endif

ifeq ($(COMP),osx_omp)
CXX = g++-5
CFLAGS = -Wall -Wextra -g -std=c++11 -O2 -fopenmp
endif


# define any directories containing header files other than /usr/include
#
INCLUDES = 

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib I'd specify
#   their path using -Lpath, something like:
LFLAGS = 

# define any libraries to link into executable:
#   if I want to link in libraries (libx.so or libx.a) I use the -llibname 
#   option, something like (this will link in libmylib.so and libm.so:
LIBS = -lo2scl_hdf -lo2scl -lhdf5

# define the C source files
SRCS = ns.cpp geo.cpp cppe.cpp main.cpp

# define the C object files 
#
# This uses Suffix Replacement within a macro:
#   $(name:string1=string2)
#         For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCS:.cpp=.o)

# define the executable file 
MAIN = cppe

#
# The following part of the makefile is generic; it can be used to 
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean

all:    $(MAIN)
	@echo  Compiler named cppe has been created.

$(MAIN): $(OBJS) 
	$(CXX) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file) 
# (see the gnu make manual section about automatic variables)
.cpp.o:
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

test:
	./cppe -run ns.in -init run0

# DO NOT DELETE THIS LINE -- make depend needs it
