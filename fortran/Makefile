# Set up conditional compilation so that one Makefile can handle
# multiple platforms/compilers. Branch according to an environmental
# variable F90. I wish someone would show me a better way of doing this.

ifeq (${F90},gfortran)  # free fortran compiler
  FFLAGS =  -g -fbounds-check -fPIC -pedantic -Wall -ffree-line-length-none
  FOUND = true
endif

ifeq (${F90},ifort)  # Intel compiler
  ifeq (${DEBUG},false)
     FFLAGS =  -O3 -fPIC
     FOUND = true
  endif
  ifeq (${DEBUG},true)
     FFLAGS =  -g -error-limit 7 -traceback -check bounds -check uninit -check pointers -fp-stack-check  -warn -heap-arrays -fpic
     FOUND = true
  endif
  ifeq (${DEBUG},)
     FFLAGS =  -g -error-limit 7 -traceback -check bounds -warn -heap-arrays -fpic
     FOUND = true
  endif

endif

ifeq (${F90},ifc)  # Intel compiler
  F90 =  /opt/intel/fc/9.1.040/bin/ifort
  FFLAGS =  -g -fPIC -error-limit 7 -traceback -check bounds -warn -heap-arrays
  FOUND = true
endif
ifeq (${F90},xlf90) # IBM compiler
  FFLAGS = -g -C -fPIC -qsuffix=f=f90  
  FOUND = true
endif
ifeq (${F90},f95) # Absoft PPC compiler
  FFLAGS = -Rb -Rc  -z2 -ea -fPIC
# B80 entry into each subroutine  ; Rb bounds; Rc array conformance
  FOUND = true
endif
ifeq (${F90},)  # If the variable isn't set, make sure compilation bombs
error:
	echo Error: environmental variable "F90" not set! Use make F90=[gfortran|ifort]
endif
ifneq (${FOUND},true) # If it's set but no matching flags then bomb
error:	
	echo Error: makefile doesn\'t have flags for this compiler; use make F90=[gfortran|ifort]
endif

SRC = itertools.f90 classes.f90
OBJS = ${SRC:.f90=.o}


.SUFFIXES :  
.SUFFIXES :  .f .f90 .f95 .o

all: libpolya.a 

libpolya.a: itertools.o classes.o 
	ar ru $@ $?
	ranlib  $@

libpolya.so:  itertools.o classes.o 
	ifort -dynamiclib -fpic -o libpolya.so itertools.f90 classes.f90
#	ar ru $@ $?
#	ranlib  $@

.f95.o : 
	${F90} ${FFLAGS} -c $<

.f90.o : 
	${F90} ${FFLAGS} -c $<

CLEAN  = *.o *.mod *.a
clean : 
	rm -f ${CLEAN}

remake:
	make clean
	make

clobber: 
	rm -f  ${CLEAN}  *~ \#*
	make
