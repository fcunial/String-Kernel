OBJS = SLT.c SLT_MAWs.c SLT_single_string.c SLT_MAWs_single_string.c dbwt_queue.c indexed_DNA5_seq.c DNA5_tables.c  dbwt.c dbwt_utils.c mt19937ar.c DNA5_Basic_BWT.c sais.c ../malloc_count-master/malloc_count.c ../malloc_count-master/stack_count.c naive_MAWs.c 

HDRS = SLT.h dbwt_queue.h indexed_DNA5_seq.h dbwt.h dbwt_utils.h mt19937ar.h DNA5_Basic_BWT.h  SLT_MAWs.h ../malloc_count-master/malloc_count.h ../malloc_count-master/stack_count.h naive_MAWs.h 





LIBS = -ldl -lm
#CFLAGS =  -g -Wall -O2 -fopenmp
CFLAGS =  -Wall -O3 -fopenmp

CC = gcc

test_SLT_MAWs : $(OBJS) test_SLT_MAWs.o ../malloc_count-master/malloc_count.o ../malloc_count-master/stack_count.o $(HDRS) 
	$(CC) $(CFLAGS) $(OBJS) test_SLT_MAWs.o -o test_SLT_MAWs $(LIBS)
test_SLT_MAWs.o: test_SLT_MAWs.c $(HDRS) 
	$(CC) $(CFLAGS) -c test_SLT_MAWs.c



#other targets
clean:
	rm *.o

