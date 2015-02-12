# makefile directives: make -j4 COMPILER=intel THREAD=yes PAPI=yes IBMHPM=yes OPT=3 USEMETIS=yes USEGREEDY=yes
# examples  : make -j4 COMPILER=gfortran
#           : make -j4 OPT=2 PAPI=no
# -------------------
# -j4 specifies that 4 processes are available for parallel make 
# You can select a compiler by using make COMPILER=intel or COMPILER=gfortran or defining a machine name below
COMPILER=$(shell  \
  if test `dnsdomainname | grep ne.anl`; then echo intel; fi; \
  if test `dnsdomainname | grep ftd`;    then echo bgq; fi; \
  if test `dnsdomainname | grep fst`;    then echo bgq; fi; \
  )
# You can enable or disable OPENMP threading by using "THREAD=no"
# You can enable or disable the two PAPI options by using PAPIFLOP and PAPISET1. They cannot both be used
# You can change the optimization settings by using "OPT=2"
THREAD=yes
PAPI=no
IBMHPM=no
OPT=3
DEBUG=
USEMETIS=no
USEGREEDY=no
METISLIB=

# Initialization of null variables
FF=echo
LD=echo
PAPIFLAGS=
PAPIINC=
PAPILIB=
IBMFLAGS=
IBMLIB=
# ------------------------------
# INTEL INTEL INTEL INTEL
# ------------------------------
ifeq (${COMPILER}, intel)
   FF     = ifort
   LD     = ${FF}
   LDFLAGS = -O${OPT} -Bstatic ${DEBUG}
ifeq (${THREAD},yes)
   LDFLAGS = -O${OPT} -openmp -DWITHOMP -Bstatic  ${DEBUG}
endif
ifeq (${PAPI},yes)
   PAPIINC= -I/software/papi-5.1.1/O-intel-12.1/include
   PAPILIB= /software/papi-5.1.1/O-intel-12.1/lib/libpapi.a  /software/papi-5.1.1/O-intel-12.1/lib/libpfm.a
endif
endif
# ------------------------------
# GFORTRAN GFORTRAN GFORTRAN 
# ------------------------------
ifeq (${COMPILER}, gfortran)
   FF     = gfortran
   LD     = ${FF}
   LDFLAGS = -O${OPT} 
ifeq (${THREAD},yes)
   LDFLAGS = -O${OPT} -fopenmp -DWITHOMP
endif
endif
ifeq (${USEGREEDY},yes)
   EXTFLAGS=-DUSEGREEDY
   METISLIB=
endif
ifeq (${USEMETIS},yes)
   EXTFLAGS=-DUSEMETIS
   METISLIB=../metis/libmetis.a
endif
# ------------------------------
# IBM IBM IBM IBM IBM IBM IBM
# ------------------------------
ifeq (${COMPILER}, bgq)
   FF     = mpixlf90_r
   LD     = ${FF}
   LDFLAGS = -O${OPT} ${DEBUG}
ifeq (${THREAD},yes)
   LDFLAGS = -O${OPT} -qsmp=omp -WF,-DWITHOMP ${DEBUG}
endif
ifeq (${IBMHPM},yes)
   IBMFLAGS= -WF,-DWITHBGQHPM
# Installed version appears to have issues as the kernel counters keep getting reassigned
#   IBMLIB = -L/soft/perftools/hpctw/lib -lmpihpm  -L/bgsys/drivers/ppcfloor/bgpm/lib -lbgpm 
# Use Vitali's stuff
   IBMLIB = -L//home/morozov/lib -lmpihpm  -L/bgsys/drivers/ppcfloor/bgpm/lib -lbgpm 
endif
ifeq (${USEGREEDY},yes)
   EXTFLAGS=-WF,-DUSEGREEDY
   METISLIB=
endif
ifeq (${USEMETIS},yes)
   EXTFLAGS=-WF,-DUSEMETIS
   METISLIB=../metis/libmetis.a
endif
endif


# ------------------------------
# Add the PAPI preprocessing stuff
TARGET = SNaCFE.x
ifeq (${PAPI},yes)
   PAPIFLAGS = -DWITHPAPI ${PAPIINC}
   TARGET = SNaCFE_PAPI.x
endif
FFLAGS = -c ${LDFLAGS} 

OBJS =  Basic_CommonBlock.o Basic_Abort.o Basic_CheckError.o \
        CommonBlock.o Method_Krylov.o GenerateXb.o PrintSummary.o Verify.o Zero.o \
        Import_PMesh.o BuildAngleCubature.o AssembleNZmatrix.o StencilNZmatrix.o \
        SolveWGS.o FGMRES_Threaded.o \
        SolveWGS_PassThrough_AVE1.o \
        SolveWGS_PassThrough_AVE2.o \
        SolveWGS_PassThrough_AVE3.o \
        SolveWGS_PassThrough_PC.o \
        ApplyA_AVE1_Tet_SUPG.o ApplyA_AVE2_Tet_SUPG.o
all: $(TARGET)

print_vars: 
	echo "host: " ${COMPILER}
	echo "comp: " ${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS}
	echo "link: " ${LD} ${LDFLAGS} ${PAPIFLAGS}
	echo "papi flag : " ${PAPI} 

CommonBlock.o: CommonBlock.F90 PROTEUS_Preprocess.h
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o CommonBlock.o CommonBlock.F90
Method_Krylov.o: Method_Krylov.F90
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o Method_Krylov.o Method_Krylov.F90
Basic_CommonBlock.o: Basic_CommonBlock.F90
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o Basic_CommonBlock.o Basic_CommonBlock.F90

Basic_Abort.o: Basic_Abort.F90 Basic_CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o Basic_Abort.o Basic_Abort.F90
Basic_CheckError.o: Basic_CheckError.F90 Basic_CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o Basic_CheckError.o Basic_CheckError.F90

GenerateXb.o: GenerateXb.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o GenerateXb.o GenerateXb.F90
PrintSummary.o: PrintSummary.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o PrintSummary.o PrintSummary.F90
Verify.o: Verify.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o Verify.o Verify.F90
Zero.o: Zero.F90 
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o Zero.o Zero.F90
Import_PMesh.o: Import_PMesh.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o Import_PMesh.o Import_PMesh.F90

BuildAngleCubature.o: BuildAngleCubature.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o BuildAngleCubature.o BuildAngleCubature.F90
StencilNZmatrix.o: StencilNZmatrix.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o StencilNZmatrix.o StencilNZmatrix.F90
AssembleNZmatrix.o: AssembleNZmatrix.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o AssembleNZmatrix.o AssembleNZmatrix.F90

SolveWGS.o: SolveWGS.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o SolveWGS.o SolveWGS.F90
Method_FGMRES.o: Method_FGMRES.F90 CommonBlock.o Method_Krylov.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o Method_FGMRES.o Method_FGMRES.F90
FGMRES_Threaded.o: FGMRES_Threaded.F90 CommonBlock.o Method_Krylov.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o FGMRES_Threaded.o FGMRES_Threaded.F90

SolveWGS_PassThrough_AVE1.o: SolveWGS_PassThrough_AVE1.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o SolveWGS_PassThrough_AVE1.o SolveWGS_PassThrough_AVE1.F90
SolveWGS_PassThrough_AVE2.o: SolveWGS_PassThrough_AVE2.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o SolveWGS_PassThrough_AVE2.o SolveWGS_PassThrough_AVE2.F90
SolveWGS_PassThrough_AVE3.o: SolveWGS_PassThrough_AVE3.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o SolveWGS_PassThrough_AVE3.o SolveWGS_PassThrough_AVE3.F90
SolveWGS_PassThrough_PC.o: SolveWGS_PassThrough_PC.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o SolveWGS_PassThrough_PC.o SolveWGS_PassThrough_PC.F90

ApplyA_AVE1_Tet_SUPG.o: ApplyA_AVE1_Tet_SUPG.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o ApplyA_AVE1_Tet_SUPG.o ApplyA_AVE1_Tet_SUPG.F90
ApplyA_AVE2_Tet_SUPG.o: ApplyA_AVE2_Tet_SUPG.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o ApplyA_AVE2_Tet_SUPG.o ApplyA_AVE2_Tet_SUPG.F90

${TARGET}: ${OBJS} Driver.F90
	$(LD) ${LDFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o ${TARGET} ${OBJS} Driver.F90 ${IBMLIB} ${PAPILIB} ${METISLIB}

clean:
	rm -f ${OBJS} *.mod
