# makefile directives: make -j4 COMPILER=intel
# examples  : make -j4 COMPILER=gfortran
#           : make -j4 COMPILER=intel
# -------------------
# -j4 specifies that 4 processes are available for parallel make 
# You can select a compiler by using make COMPILER=intel or COMPILER=gfortran or defining a machine name below
COMPILER=intel
COMPILER=$(shell  \
  if test `dnsdomainname | grep ne.anl`; then echo intel; fi; \
  )

# Initialization of null variables
FF=echo
LD=echo
# ------------------------------
# INTEL INTEL INTEL INTEL
# ------------------------------
ifeq (${COMPILER}, intel)
   FF     = ifort
   LD     = ${FF}
   LDFLAGS = 
endif
# ------------------------------
# GFORTRAN GFORTRAN GFORTRAN 
# ------------------------------
ifeq (${COMPILER}, gfortran)
   FF     = gfortran
   LD     = ${FF}
   LDFLAGS = 
endif

# ------------------------------
TARGET = makemesh.x
FFLAGS = -c ${LDFLAGS} 

OBJS =  CommonBlock.o Build_Grid_Tet.o Export_Mesh.o Transform_Hex20.o
all: $(TARGET)

CommonBlock.o: CommonBlock.F90 PROTEUS_Preprocess.h
	${FF} ${FFLAGS} -o CommonBlock.o CommonBlock.F90
Build_Grid_Tet.o: Build_Grid_Tet.F90 CommonBlock.o
	${FF} ${FFLAGS} -o Build_Grid_Tet.o Build_Grid_Tet.F90
Export_Mesh.o: Export_Mesh.F90 CommonBlock.o
	${FF} ${FFLAGS} -o Export_Mesh.o Export_Mesh.F90
Transform_Hex20.o: Transform_Hex20.F90 
	${FF} ${FFLAGS} -o Transform_Hex20.o Transform_Hex20.F90

${TARGET}: ${OBJS} Driver.F90
	$(LD) ${LDFLAGS} -o ${TARGET} ${OBJS} Driver.F90

clean:
	rm -f ${OBJS} *.mod
