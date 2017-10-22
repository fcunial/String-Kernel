CC=gcc
CFLAGS=-Wall -O3 -fopenmp
LIBS=-ldl -lm


all: $(PROGRAM_1)




# ---- MAIN PROGRAMS ----

PROGRAM_1=run_MAWs_single.c
program-1: dbwt malloc-count random iterator maws_single io $(PROGRAM_1)
	$(CC) $(CFLAGS) $(PROGRAM_1) $(LIBS) -o run_MAWs_single  




# ---- COMPONENTS ----

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


IO_DIR=./io
IO_OBJS=$(IO_DIR)/io.c
IO_HDRS=$(IO_DIR)/io.h
io: $(IO_OBJS) $(IO_HDRS)
	$(CC) $(CFLAGS) -c $(IO_OBJS) $(IO_HDRS)




# ---- CALLBACKS ----

CALLBACKS_DIR=./callbacks

MAWS_SINGLE_OBJS=$(CALLBACKS_DIR)/MAWs_single.c $(MAWS_SINGLE_DIR)/io.c $(MAWS_SINGLE_DIR)/run_MAWs_single.c
MAWS_SINGLE_HDRS=$(CALLBACKS_DIR)/MAWs_single.h
maws_single: $(MAWS_SINGLE_OBJS) $(MAWS_SINGLE_HDRS)
	$(CC) $(CFLAGS) -c $(MAWS_SINGLE_OBJS) $(MAWS_SINGLE_HDRS)




# ---- CLEANING ----

clean:
	rm $(CALLBACKS_DIR)/*.o $(IO_DIR)/*.o $(ITERATOR_DIR)/*.o $(RANDOM_DIR)/*.o $(MALLOC_COUNT_DIR)/*.o $(DBWT_DIR)/*.o
 
