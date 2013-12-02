#.SILENT:

.SUFFIXES:   .o .f 

FFLAGS = -ffixed-line-length-none
F77      = gfortran
LINK	 = gfortran
LINKLIB  =  -lm -lgfortran
CERNLINK = -L/cern/cernlib/2006b/x86_64/lib -lmathlib -lkernlib -lpacklib


# ==== make-rules =========================================================


.f.o :
	echo 'Compiling $<'
	$(F77) $(FFLAGS) -c $<


# ==== executable programs ===============================================
all:	ugcalv2ud

ugcalv2ud : ugcalv2ud.o
		echo Linking ugcalv2ud
		$(LINK) ugcalv2ud.o $(CERNLINK) $(LINKLIB) -o ugcalv2ud

clean :
	rm -f *.o ugcalv2ud
