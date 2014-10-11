SHELL=/bin/csh
.SUFFIXES: .dia .tst .vl

SCIDIR=../../../

TESTS= linprog.dia examples.dia

VL=$(ALG:.dia=.vl) $(TESTS:.dia=.vl) 

all : algebra basic 

vl  : $(VL)

algebra : $(ALG) 

basic : $(TESTS) 

tests: all

.tst.dia :
	@echo -n "testing $*: "
	@$(SCIDIR)/bin/nsp -nw -e "nsptest('$*.tst');quit" 

#export G_DEBUG=gc-friendly;export G_SLICE=always-malloc;valgrind --suppressions=tests/valgrind.supp --leak-check=full  bin/scilex -nw 

.tst.vl :
	@echo -n "valgrind $*: "
	@setenv G_DEBUG gc-friendly; setenv G_SLICE always-malloc;valgrind --suppressions=valgrind.supp --leak-check=full  ../bin/scilex -nw -e "exec('$*.tst');quit" >& $*.vl
	@echo  " $*.vl"

distclean	:
	@echo cleaning tests
	@$(RM) $(TESTS) $(ALG) $(VL) 
	@$(SCIDIR)/bin/nsplibtool OBJS=addinter.o LIBRARY=libaddinter distclean

clean:	distclean

