SHELL = /bin/sh

SCIDIR=../../..
SCIDIR1=..\..\..

LIBRARY=libcplex.a

CPLEX=/opt/ibm/ILOG/CPLEX_Studio126/cplex
CPLEX_INC= -I /opt/ibm/ILOG/CPLEX_Studio126/cplex/include 
CPLEX_LIB= -L/opt/ibm/ILOG/CPLEX_Studio126/cplex/bin/x86_linux/ -Wl,-R$(CPLEX)/bin/x86_linux -lcplex1260 -lpthread 

OBJS= nspcplex-In.o nspcplex.o stochdeclp.o stochdeclp_cplex.o

include $(SCIDIR)/Makefile.incl

CFLAGS = $(CC_OPTIONS) $(CPLEX_INC)

# extra libraries needed for linking 
# it is mandatory on win32 to give this extra argument.
OTHERLIBS=$(CPLEX_LIB)

include $(SCIDIR)/config/Makeso.incl

all:: Makefile.mak

Makefile.mak	: Makefile
	$(SCIDIR)/scripts/Mak2VCMak Makefile

Makefile.amk	: Makefile
	$(SCIDIR)/scripts/Mak2ABSMak Makefile

distclean:: clean 

clean:: 
	@$(RM) *.o *.lo 
	@$(RM) -r */.libs 

