LOC_GPP   = /usr/bin/g++
MPI_LIBS  = -lm -lpthread -lfftw3
MPI_OBJS  = mpi_ecosim_threads.o mpi_csvparser.o mpi_gauss_threads.o mpi_nrutil.o
MPI_MAIN  = ecosim_threads.cpp ecosim_threads.h ecosim_UserInputs.h \
            simio.cpp fitfunctions.cpp sense.cpp simlabs.cpp fftest.cpp MPI_sim.cpp 
MPI_FLAGS = -O3 -I/home/aydink/include -L/home/aydink/lib -D MPI_COMPILE

################################################################################

LOC_LIBS  = -lm -lgsl -lgslcblas -lpthread -lfftw3 
#LOC_LIBS  = -lm -lgsl -lgslcblas -lpthread 
LOC_OBJS  = loc_ecosim_threads.o loc_csvparser.o loc_gauss_threads.o loc_nrutil.o
LOC_MAIN  = ecosim_threads.cpp ecosim_threads.h ecosim_UserInputs.h \
            simio.cpp fitfunctions.cpp sense.cpp simlabs.cpp fftest.cpp 
#LOC_FLAGS = -O3 -static  -ggdb
LOC_FLAGS = -O3 -g

################################################################################

all :
	@echo "Use 'make mpi' or 'make no-mpi'"

################################################################################  
mpi : ${MPI_OBJS}
	@echo "Making mpi version"
	/opt/openmpi/bin/mpiCC -o ecosim ${MPI_FLAGS} ${MPI_OBJS} ${MPI_LIBS}
  
mpi_ecosim_threads.o : ${MPI_MAIN} 
	/opt/openmpi/bin/mpiCC -c ${MPI_FLAGS} ecosim_threads.cpp
	mv ecosim_threads.o mpi_ecosim_threads.o
	
mpi_csvparser.o: csvparser.cpp csvparser.h
	/opt/openmpi/bin/mpiCC -c ${MPI_FLAGS} csvparser.cpp
	mv csvparser.o mpi_csvparser.o
	
mpi_gauss_threads.o: gauss_threads.cpp gauss_threads.h RandomSeeds.h
	/opt/openmpi/bin/mpiCC -c ${MPI_FLAGS} gauss_threads.cpp
	mv gauss_threads.o mpi_gauss_threads.o
	
mpi_nrutil.o: nrutil.cpp nrutil.h
	/opt/openmpi/bin/mpiCC -c ${MPI_FLAGS} nrutil.cpp
	mv nrutil.o mpi_nrutil.o 

################################################################################
  
no-mpi : ${LOC_OBJS}
	@echo "Making no-mpi verion"
	${LOC_GPP} -o ecosim ${LOC_FLAGS} ${LOC_OBJS} ${LOC_LIBS}
#	rm ${OBJS}
  
loc_ecosim_threads.o : ${LOC_MAIN} 
	${LOC_GPP} -c ${LOC_FLAGS} ecosim_threads.cpp
	mv ecosim_threads.o loc_ecosim_threads.o

loc_csvparser.o: csvparser.cpp csvparser.h
	${LOC_GPP} -c ${LOC_FLAGS} csvparser.cpp
	mv csvparser.o loc_csvparser.o
	
loc_gauss_threads.o: gauss_threads.cpp gauss_threads.h RandomSeeds.h
	${LOC_GPP} -c ${LOC_FLAGS} gauss_threads.cpp
	mv gauss_threads.o loc_gauss_threads.o 
	
loc_nrutil.o: nrutil.cpp nrutil.h
	${LOC_GPP} -c ${LOC_FLAGS} nrutil.cpp
	mv nrutil.o loc_nrutil.o
