SHELL = /bin/sh

SCIDIR=../../..
SCIDIR1=..\..\..

LIBRARY=libcplex.a

CPLEX=/opt/ibm/ILOG/CPLEX_Studio201/cplex
CPLEXLIB=-lcplex2010

CPLEX_INC= -I $(CPLEX)/include 
CPLEX_LIB= -L$(CPLEX)/bin/x86-64_linux/ -Wl,-R$(CPLEX)/bin/x86-64_linux $(CPLEXLIB)

OBJS= nspcplex-In.o nspcplex.o 

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

