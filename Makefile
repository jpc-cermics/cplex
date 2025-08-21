SHELL = /bin/sh

SILENT=$(findstring -s,$(MFLAGS))
NSP=../../bin/nsp 
GNUMAKEFLAGS=--no-print-dir

all:
	@if test -f Path.incl; then \
		$(MAKE) $(MFLAGS) all-dirs ; \
	else \
	  ( if test -f $(NSP); then \
		( if test "x$(SILENT)" != "x-s"; then echo "running builder"; fi && \
		$(NSP) -nw -ns -e "exec('builder.sce');quit" -errcatch > /dev/null ) ; \
	    else \
	      echo "Fisrt time you run make;start nsp and run the file builder.sce"; \
	    fi);\
	fi

SUBDIRS =src
DIR=

all-dirs:
	@case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS) ;\
	do \
		(cd $$i && if test "x$(SILENT)" != "x-s"; then echo "making all in $(DIR)$$i ";fi && \
		$(MAKE) $(MFLAGS) DIR=$(DIR)$$i/ all ); \
	   	IER=$$? &&\
	   	case $$IER in\
	    	0) ;;\
	    	*) echo "make $@ in sub directory $$i failed"; \
	       	   case '${MFLAGS}' in *[k]*) echo "carrying on compilation (-k used)";; *) exit $$IER;;esac;\
	   	esac;\
	done

clean distclean ::
	@if test -f Path.incl; then \
	(case '${MFLAGS}' in *[ik]*) set +e;; esac; \
	for i in $(SUBDIRS) ;\
	do \
		(cd $$i && if test "x$(SILENT)" != "x-s"; then echo "making $@ in $(DIR)$$i ";fi && \
		$(MAKE) $(MFLAGS) DIR=$(DIR)$$i/ $@ ); \
	   	IER=$$? &&\
	   	case $$IER in\
	    	0) ;;\
	    	*) echo "make $@ in sub directory $$i failed"; \
	       	   case '${MFLAGS}' in *[k]*) echo "carrying on compilation (-k used)";; *) exit $$IER;;esac;\
	   	esac;\
	done); \
	else $(MAKE) $(MFLAGS) distclean-base; \
	fi

distclean::
	@$(RM) Path.incl

# target to make a distclean when we do not have a Path.incl

distclean-base::
	@find . \( -name .libs -o -name '*.o' -o -name '*.so' -o -name '*.a' -name '*.bin' \) \
		-exec \rm -f {} \;

