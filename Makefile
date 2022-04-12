SHELL=/bin/sh
MFLAGS=-s --no-print-directory

include Path.incl
include $(SCIDIR)/Makefile.incl

all :: builder.sce
	@echo "running builder (be patient)"
	@$(SCIDIR)/bin/nsp -nw -e "exec('builder.sce');quit" -errcatch > /dev/null 2>&1;
	@echo "At prompt, enter:";
	@echo "-->exec loader.sce";
	@echo "----------------------------------------------------";

all ::
	cd src && $(MAKE) all

clean ::
	@echo "Clean src "
	@cd src && $(MAKE) $(MFLAGS) clean

cleanm ::
	@echo "Clean macros"
	@cd macros && $(MAKE) $(MFLAGS) clean

distclean ::
	@echo "Clean src"
	@cd src && $(MAKE) distclean > /dev/null 2>&1;

cleanm	::
	@echo "Clean macros"
	@cd macros && $(MAKE) distclean > /dev/null 2>&1;

PATH_INCL= $(wildcard Path.incl)

test ::
ifeq ($(PATH_INCL),Path.incl)
	@echo $(PATH_INCL) "is already present"
else
	@echo Path.incl "copied from ../../"
	@cp ../../Path.incl .
endif

up:
	@echo "Update cvs"
	@cvs up 

