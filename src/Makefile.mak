#-----------------------------
# generated from Makefile: DO NOT EDIT
# -----------------------------
SHELL = /bin/sh

LIBRARY=libcplex.lib

CPLEX=/opt/ibm/ILOG/CPLEX_Studio201/cplex
CPLEXLIB=-lcplex2010
CPLEX=/opt/ibm/ILOG/CPLEX_Studio221/cplex
CPLEXLIB=-lcplex2210

CPLEX_INC= -I $(CPLEX)/include 
CPLEX_LIB= -L$(CPLEX)/bin/x86-64_linux/ -Wl,-R$(CPLEX)/bin/x86-64_linux $(CPLEXLIB)

OBJS= nspcplex-In.obj nspcplex.obj 

include ../Path.incl
include $(SCIDIR)/Makefile.incl.mak

CFLAGS = $(CC_OPTIONS) $(CPLEX_INC)

# extra libraries needed for linking 
# it is mandatory on win32 to give this extra argument.
OTHERLIBS=$(CPLEX_LIB)

include $(SCIDIR)/config/Makeso.incl



Makefile.mak	: Makefile
	$(SCIDIR)/scripts/Mak2VCMak Makefile

Makefile.libmk	: Makefile
	$(SCIDIR)/scripts/Mak2ABSMak Makefile

distclean:: clean 

clean:: 
	@$(RM) *.obj *.lo 
	@$(RM) -r */.libs 

