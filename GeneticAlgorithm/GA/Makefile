COMP = gcc
LIB = -lgsl -lgslcblas -lm
OPT = -Wall -O3 -DHAVE_INLINE
OBJECTS = main.o find_eigen.o read_write_matrix.o random_number.o mutate.o

FindArrangGA: ${OBJECTS}
	${COMP} ${OBJECTS} -o FindArrangGA ${OPT} ${LIB}

read_write_matrix.o: common.h read_write_matrix.c 
	${COMP} read_write_matrix.c -c ${OPT}

random_number.o: common.h random_number.c 
	${COMP} random_number.c -c ${OPT}

find_eigen.o: common.h find_eigen.c 
	${COMP} find_eigen.c -c ${OPT}

mutate.o: common.h mutate.c 
	${COMP} mutate.c -c ${OPT}


main.o: main.c common.h find_eigen.h read_write_matrix.h random_number.h mutate.h
	${COMP} main.c -c ${OPT}

clean: 
	rm -rf *.o *~ *.log 
