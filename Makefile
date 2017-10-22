




DBWT_DIR=./dbwt
DBWT_OBJS=$(DBWT_DIR)/dbwt.c $(DBWT_DIR)/dbwt_queue.c $(DBWT_DIR)/dbwt_utils.c
DBWT_HDRS=$(DBWT_HDRS)/dbwt.h $(DBWT_HDRS)/dbwt_queue.h $(DBWT_HDRS)/dbwt_utils.h
dbwt: $(DBWT_OBJS) $(DBWT_HDRS)
	$(CC) $(CFLAGS) -c $(DBWT_OBJS) $(LIBS) 

MALLOC_COUNT_DIR=./malloc_count
MALLOC_COUNT_OBJS=$(MALLOC_COUNT_DIR)/malloc_count.c $(MALLOC_COUNT_DIR)/stack_count.c 
MALLOC_COUNT_HDRS=$(MALLOC_COUNT_DIR)/malloc_count.h $(MALLOC_COUNT_DIR)/stack_count.h
malloc-count: $(MALLOC_COUNT_OBJS) $(MALLOC_COUNT_HDRS)
	$(CC) $(CFLAGS) -c $(MALLOC_COUNT_OBJS) $(LIBS)

RANDOM_DIR=./random
RANDOM_OBJS=$(RANDOM_DIR)/mt19937ar.c
RANDOM_HDRS=$(RANDOM_DIR)/mt19937ar.h
random: $(RANDOM_OBJS) $(RANDOM_HDRS)
	$(CC) $(CFLAGS) -c $(RANDOM_OBJS) $(LIBS)

ITERATOR_DIR=./iterator
ITERATOR_OBJS=$(ITERATOR_DIR)/DNA5_tables.c $(ITERATOR_DIR)/indexed_DNA5_seq.c $(ITERATOR_DIR)/DNA5_Basic_BWT.c $(ITERATOR_DIR)/SLT_single_string.c 
ITERATOR_HDRS=$(ITERATOR_DIR)/indexed_DNA5_seq.h $(ITERATOR_DIR)/DNA5_Basic_BWT.h $(ITERATOR_DIR)/SLT_single_string.h
iterator: $(ITERATOR_OBJS) $(ITERATOR_HDRS)
	$(CC) $(CFLAGS) -c $(ITERATOR_OBJS) $(LIBS)

MAWS_SINGLE_DIR=./ParallelGeneralized
MAWS_SINGLE_OBJS=$(MAWS_SINGLE_DIR)/MAWs_single.c $(MAWS_SINGLE_DIR)/io.c $(MAWS_SINGLE_DIR)/run_MAWs_single.c



SLT_MAWs.c  SLT_MAWs_single_string.c mt19937ar.c  sais.c  naive_MAWs.c 

HDRS = SLT.h  mt19937ar.h  SLT_MAWs.h  naive_MAWs.h 






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

