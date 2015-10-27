#   @(#)Makefile	4.1 7/9/93 
#
#  Makefile to assemble MesaLib in the directory mesa/library
#  Machine dependent portion of the library.
#
OPSYS = sun
#
all: $(OPSYS)
.RECURSIVE: $(OPSYS)
$(OPSYS): FORCE
	cd $@; $(MAKE) $(MFLAGS)

sources: FORCE
	cd unicos ;    $(MAKE) $(MFLAGS) source
	cd aix ;       $(MAKE) $(MFLAGS) source
	cd stardent ;  $(MAKE) $(MFLAGS) source
	cd sun ;       $(MAKE) $(MFLAGS) source
	cd ultrix ;    $(MAKE) $(MFLAGS) source
	cd hp ;        $(MAKE) $(MFLAGS) source
	cd sgi;        $(MAKE) $(MFLAGS) source

clean: FORCE
	cd unicos ;    $(MAKE) $(MFLAGS) clean
	cd aix ;       $(MAKE) $(MFLAGS) clean
	cd stardent ;  $(MAKE) $(MFLAGS) clean
	cd sun ;       $(MAKE) $(MFLAGS) clean
	cd ultrix ;    $(MAKE) $(MFLAGS) clean
	cd hp;	       $(MAKE) $(MFLAGS) clean
	cd sgi;	       $(MAKE) $(MFLAGS) clean

link: FORCE
	cd unicos ;    $(MAKE) $(MFLAGS) link
	cd aix ;       $(MAKE) $(MFLAGS) link
	cd stardent ;  $(MAKE) $(MFLAGS) link
	cd sun ;       $(MAKE) $(MFLAGS) link
	cd ultrix ;    $(MAKE) $(MFLAGS) link
	cd hp ;        $(MAKE) $(MFLAGS) link
	cd sgi;        $(MAKE) $(MFLAGS) link

FORCE:
