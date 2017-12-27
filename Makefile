CC=/usr/bin/gcc
CFLAGS=-Wall -O3 #-fopenmp
LIBS=-ldl -lm
ROOT_DIR=/Users/ramseysnow/git/String-Kernel
.PHONY: all clean   program-1   dbwt malloc-count random iterator io maws-single


all: program-1
	




# ---- MAIN PROGRAMS ----

PROGRAMS=$(PROGRAM_1)

PROGRAM_1=$(ROOT_DIR)/run_MAWs_single
program-1: $(PROGRAM_1).c io dbwt malloc-count iterator maws-single
	$(CC) $(CFLAGS) $(LIBS) $(PROGRAM_1).c $(IO_OBJS) $(DBWT_OBJS) $(MALLOC_COUNT_OBJS) $(ITERATOR_OBJS) $(MAWS_SINGLE_OBJS) -o $(PROGRAM_1)




# ---- COMPONENTS ----

DBWT_DIR=$(ROOT_DIR)/dbwt
DBWT_SRC=$(DBWT_DIR)/dbwt.c $(DBWT_DIR)/dbwt_queue.c $(DBWT_DIR)/dbwt_utils.c $(DBWT_DIR)/sais.c
DBWT_HDRS=$(DBWT_DIR)/dbwt.h $(DBWT_DIR)/dbwt_queue.h $(DBWT_DIR)/dbwt_utils.h
DBWT_OBJS=$(DBWT_DIR)/dbwt.o $(DBWT_DIR)/dbwt_queue.o $(DBWT_DIR)/dbwt_utils.o $(DBWT_DIR)/sais.o
dbwt: $(DBWT_SRC) $(DBWT_HDRS)
	cd $(DBWT_DIR) && $(CC) $(CFLAGS) -c *.c


MALLOC_COUNT_DIR=$(ROOT_DIR)/malloc_count
MALLOC_COUNT_SRC=$(MALLOC_COUNT_DIR)/malloc_count.c $(MALLOC_COUNT_DIR)/stack_count.c 
MALLOC_COUNT_HDRS=$(MALLOC_COUNT_DIR)/malloc_count.h $(MALLOC_COUNT_DIR)/stack_count.h
MALLOC_COUNT_OBJS=$(MALLOC_COUNT_DIR)/malloc_count.o $(MALLOC_COUNT_DIR)/stack_count.o 
malloc-count: $(MALLOC_COUNT_SRC) $(MALLOC_COUNT_HDRS)
	cd $(MALLOC_COUNT_DIR) && $(CC) $(CFLAGS) -c *.c


RANDOM_DIR=$(ROOT_DIR)/random
RANDOM_SRC=$(RANDOM_DIR)/mt19937ar.c
RANDOM_HDRS=$(RANDOM_DIR)/mt19937ar.h
RANDOM_OBJS=$(RANDOM_DIR)/mt19937ar.o
random: $(RANDOM_SRC) $(RANDOM_HDRS)
	cd $(RANDOM_DIR) && $(CC) $(CFLAGS) $(LIBS) -c $(RANDOM_SRC)


ITERATOR_DIR=$(ROOT_DIR)/iterator
ITERATOR_SRC=$(ITERATOR_DIR)/DNA5_tables.c $(ITERATOR_DIR)/indexed_DNA5_seq.c $(ITERATOR_DIR)/DNA5_Basic_BWT.c $(ITERATOR_DIR)/SLT_single_string.c 
ITERATOR_HDRS=$(ITERATOR_DIR)/indexed_DNA5_seq.h $(ITERATOR_DIR)/DNA5_Basic_BWT.h $(ITERATOR_DIR)/SLT_single_string.h
ITERATOR_OBJS=$(ITERATOR_DIR)/DNA5_tables.o $(ITERATOR_DIR)/indexed_DNA5_seq.o $(ITERATOR_DIR)/DNA5_Basic_BWT.o $(ITERATOR_DIR)/SLT_single_string.o 
iterator: $(ITERATOR_SRC) $(ITERATOR_HDRS)
	cd $(ITERATOR_DIR) && $(CC) $(CFLAGS) -c $(ITERATOR_SRC)


IO_DIR=$(ROOT_DIR)/io
IO_SRC=$(IO_DIR)/io.c
IO_HDRS=$(IO_DIR)/io.h
IO_OBJS=$(IO_DIR)/io.o
io: $(IO_SRC) $(IO_HDRS)
	cd $(IO_DIR) && $(CC) $(CFLAGS) -c $(IO_SRC)




# ---- CALLBACKS ----

CALLBACKS_DIR=$(ROOT_DIR)/callbacks

MAWS_SINGLE_SRC=$(CALLBACKS_DIR)/MAWs_single.c
MAWS_SINGLE_HDRS=$(CALLBACKS_DIR)/MAWs_single.h
MAWS_SINGLE_OBJS=$(CALLBACKS_DIR)/MAWs_single.o
maws-single: $(MAWS_SINGLE_SRC) $(MAWS_SINGLE_HDRS)
	cd $(CALLBACKS_DIR) && $(CC) $(CFLAGS) -c $(MAWS_SINGLE_SRC)



# ---- CLEANING ----

clean:
	rm $(CALLBACKS_DIR)/*.o $(IO_DIR)/*.o $(ITERATOR_DIR)/*.o $(RANDOM_DIR)/*.o $(MALLOC_COUNT_DIR)/*.o $(DBWT_DIR)/*.o $(PROGRAMS)
 
