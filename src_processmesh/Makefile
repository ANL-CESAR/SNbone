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

EXTFLAGS=-DUSEMETIS
METISLIB=../metis/libmetis.a

# Initialization of null variables
FF=echo
LD=echo
OPT=3
# ------------------------------
# INTEL INTEL INTEL INTEL
# ------------------------------
ifeq (${COMPILER}, intel)
   FF     = ifort
   LD     = ${FF}
   LDFLAGS = -O${OPT} -Bstatic ${DEBUG}
endif
# ------------------------------
# GFORTRAN GFORTRAN GFORTRAN 
# ------------------------------
ifeq (${COMPILER}, gfortran)
   FF     = gfortran
   LD     = ${FF}
   LDFLAGS = -O${OPT} 
endif

TARGET = processmesh.x
FFLAGS = -c ${LDFLAGS} 

OBJS =  CommonBlock.o \
        RearrangeConn.o GetSpatialMatrices.o Element_NatTet.o Export_PMesh.o ExportVTK.o Import_Mesh.o \
        GetAdjacencyTElement.o GetVertextoThread_GREEDY.o GetVertextoThread_METIS.o \
        GetTElementOrder_GREEDY.o GetTElementOrder_METIS.o  \
        GetTVertexOrder.o CheckAdjacency.o
all: $(TARGET)

print_vars: 
	echo "host: " ${COMPILER}
	echo "comp: " ${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS}
	echo "link: " ${LD} ${LDFLAGS} ${PAPIFLAGS}
	echo "papi flag : " ${PAPI} 

CommonBlock.o: CommonBlock.F90 PROTEUS_Preprocess.h
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o CommonBlock.o CommonBlock.F90

Export_PMesh.o: Export_PMesh.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o Export_PMesh.o Export_PMesh.F90
ExportVTK.o: ExportVTK.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o ExportVTK.o ExportVTK.F90
Import_Mesh.o: Import_Mesh.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o Import_Mesh.o Import_Mesh.F90

RearrangeConn.o: RearrangeConn.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o RearrangeConn.o RearrangeConn.F90
GetSpatialMatrices.o: GetSpatialMatrices.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o GetSpatialMatrices.o GetSpatialMatrices.F90
GetAdjacencyTElement.o: GetAdjacencyTElement.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o GetAdjacencyTElement.o GetAdjacencyTElement.F90
GetVertextoThread_GREEDY.o: GetVertextoThread_GREEDY.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o GetVertextoThread_GREEDY.o GetVertextoThread_GREEDY.F90
GetVertextoThread_METIS.o: GetVertextoThread_METIS.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o GetVertextoThread_METIS.o GetVertextoThread_METIS.F90

Element_NatTet.o: Element_NatTet.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o Element_NatTet.o Element_NatTet.F90
GetTElementOrder_GREEDY.o: GetTElementOrder_GREEDY.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o GetTElementOrder_GREEDY.o GetTElementOrder_GREEDY.F90
GetTElementOrder_METIS.o: GetTElementOrder_METIS.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o GetTElementOrder_METIS.o GetTElementOrder_METIS.F90
GetTVertexOrder.o: GetTVertexOrder.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o GetTVertexOrder.o GetTVertexOrder.F90

CheckAdjacency.o: CheckAdjacency.F90 CommonBlock.o
	${FF} ${FFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o CheckAdjacency.o CheckAdjacency.F90

${TARGET}: ${OBJS} Driver.F90
	$(LD) ${LDFLAGS} ${PAPIFLAGS} ${EXTFLAGS} ${IBMFLAGS} -o ${TARGET} ${OBJS} Driver.F90 ${METISLIB}

clean:
	rm -f ${OBJS} *.mod
