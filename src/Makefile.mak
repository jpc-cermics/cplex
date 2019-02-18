#-----------------------------
# generated from Makefile: DO NOT EDIT
# -----------------------------
SHELL = /bin/sh

SCIDIR=../../..
SCIDIR1=..\..\..

LIBRARY=libcplex.lib

CPLEX=/opt/ibm/ILOG/CPLEX_Studio126/cplex
CPLEX_INC= -I /opt/ibm/ILOG/CPLEX_Studio126/cplex/include 
CPLEX_LIB= -L/opt/ibm/ILOG/CPLEX_Studio126/cplex/bin/x86-64_linux/ -Wl,-R$(CPLEX)/bin/x86-64_linux -lcplex1260 -lpthread 

OBJS= nspcplex-In.obj nspcplex.obj 

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

