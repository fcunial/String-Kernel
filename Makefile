CC="gcc"
CFLAGS=-fopenmp -Wall -O3 -m64 -mfpmath=both -march=native -mtune=native -ffast-math -ftree-vectorize -mavx -mavx2
#-fopt-info-vec-missed#-fopt-info-vec-all -g#missed - optimized - all
#-march=icelake-client -mtune=icelake-client
#-mtune=intel -mtune=generic #-masm=intel -mavx -mavx2
#-g -pg #-fno-tree-vectorize #-fopt-info-optimized-optall
#-mno-avx256-split-unaligned-load -mno-avx256-split-unaligned-store

LIB_PATH=/usr/local/lib
LIBS=-ldl -lm
DIVSUFSORT_OBJS=$(LIB_PATH)/libdivsufsort64.a
ROOT_DIR=$(CURDIR)
RUN_MAW_SINGLE=$(ROOT_DIR)/run_MAWs_single
RUN_MRW_SINGLE=$(ROOT_DIR)/run_MRWs_single
BUILDINDEX=$(ROOT_DIR)/buildIndex
BUILDBWT=$(ROOT_DIR)/buildBWT

.PHONY: all clean random


all: maws-single run_MAWs_single run_MRWs_single buildIndex buildBWT random

MALLOC_COUNT_DIR=$(ROOT_DIR)/malloc_count
MALLOC_COUNT_LIST=malloc_count.c stack_count.c malloc_count.o stack_count.o
MALLOC_COUNT_OBJS=$(MALLOC_COUNT_DIR)/malloc_count.o $(MALLOC_COUNT_DIR)/stack_count.o 

RANDOM_DIR=$(ROOT_DIR)/random
RANDOM_OBJS=$(RANDOM_DIR)/mt19937ar.o

ITERATOR_DIR=$(ROOT_DIR)/iterator
ITERATOR_LIST=DNA5_tables.c indexed_DNA5_seq.c DNA5_Basic_BWT.c SLT_single_string.c DNA5_tables.o indexed_DNA5_seq.o DNA5_Basic_BWT.o SLT_single_string.o
ITERATOR_OBJS=$(ITERATOR_DIR)/DNA5_tables.o $(ITERATOR_DIR)/indexed_DNA5_seq.o $(ITERATOR_DIR)/DNA5_Basic_BWT.o $(ITERATOR_DIR)/SLT_single_string.o

BWT_DIR=$(ROOT_DIR)/BWT
BWT_LIST=BwtBuilders.c BwtBuilders.o
BWT_OBJS=$(BWT_DIR)/BwtBuilders.o

IO_DIR=$(ROOT_DIR)/io
IO_LIST=io.c bufferedFileWriter.c bits.c io.o bufferedFileWriter.o bits.o
IO_OBJS=$(IO_DIR)/io.o $(IO_DIR)/bufferedFileWriter.o $(IO_DIR)/bits.o

CALLBACKS_DIR=$(ROOT_DIR)/callbacks
CALLBACKS_LIST = MAWs_single.c MAWs_single.o
MAWS_SINGLE_OBJS=$(CALLBACKS_DIR)/MAWs_single.o

VPATH=.:$(IO_DIR):$(MALLOC_COUNT_DIR):$(ITERATOR_DIR):$(CALLBACKS_DIR):$(RANDOM_DIR)

# ---- COMPONENTS ----

malloc_count.o: malloc_count.c malloc_count.h
	cd $(MALLOC_COUNT_DIR) && $(CC) $(CFLAGS) -c malloc_count.c
stack_count.o: stack_count.c stack_count.h
	cd $(MALLOC_COUNT_DIR) && $(CC) $(CFLAGS) -c stack_count.c

#mt19937ar.o: mt19937ar.c mt19937ar.h
#	cd $(RANDOM_DIR) && $(CC) $(CFLAGS) $(LIBS) -c $(RANDOM_SRC)

DNA5_tables.o: DNA5_tables.c
	cd $(ITERATOR_DIR) && $(CC) $(CFLAGS) -c DNA5_tables.c
indexed_DNA5_seq.o: indexed_DNA5_seq.c indexed_DNA5_seq.h
	cd $(ITERATOR_DIR) && $(CC) $(CFLAGS) -c indexed_DNA5_seq.c
DNA5_Basic_BWT.o: DNA5_Basic_BWT.c DNA5_Basic_BWT.h
	cd $(ITERATOR_DIR) && $(CC) $(CFLAGS) -c DNA5_Basic_BWT.c
SLT_single_string.o: SLT_single_string.c SLT_single_string.h
	cd $(ITERATOR_DIR) && $(CC) $(CFLAGS) -c SLT_single_string.c


BwtBuilders.o: BwtBuilders.c BwtBuilders.h
	cd $(BWT_DIR) && $(CC) $(CFLAGS) -c BwtBuilders.c

io.o: io.c io.h
	cd $(IO_DIR) && $(CC) $(CFLAGS) -c io.c
bufferedFileWriter.o: bufferedFileWriter.c bufferedFileWriter.h
	cd $(IO_DIR) && $(CC) $(CFLAGS) -c bufferedFileWriter.c
bits.o: bits.c bits.h
	cd $(IO_DIR) && $(CC) $(CFLAGS) -c bits.c

# ---- CALLBACKS ----
maws-single: MAWs_single.c MAWs_single.h
	cd $(CALLBACKS_DIR) && $(CC) $(CFLAGS) -c MAWs_single.c
MAWs_single.o: MAWs_single.c MAWs_single.h
	cd $(CALLBACKS_DIR) && $(CC) $(CFLAGS) -c MAWs_single.c
# ---- MAIN PROGRAMS ----

run_MAWs_single: Makefile scores.c $(RUN_MAW_SINGLE).c $(IO_LIST) $(MALLOC_COUNT_LIST) $(ITERATOR_LIST) $(CALLBACK_LIST)
		$(CC) $(CFLAGS) scores.c $(RUN_MAW_SINGLE).c $(IO_OBJS) $(MALLOC_COUNT_OBJS) $(ITERATOR_OBJS) $(MAWS_SINGLE_OBJS) $(LIBS) -o $(RUN_MAW_SINGLE)

run_MRWs_single: Makefile scores.c $(RUN_MRW_SINGLE).c $(IO_LIST) $(MALLOC_COUNT_LIST) $(ITERATOR_LIST) $(CALLBACK_LIST)
	$(CC) $(CFLAGS) scores.c $(RUN_MRW_SINGLE).c $(IO_OBJS) $(MALLOC_COUNT_OBJS) $(ITERATOR_OBJS) $(MAWS_SINGLE_OBJS) $(LIBS) -o $(RUN_MRW_SINGLE)

buildIndex: Makefile buildIndex.c $(BUILDINDEX).c $(IO_LIST) $(MALLOC_COUNT_LIST) $(ITERATOR_LIST)
		$(CC) $(CFLAGS) $(BUILDINDEX).c $(IO_OBJS) $(MALLOC_COUNT_OBJS) $(ITERATOR_OBJS) $(LIBS)  -o $(BUILDINDEX)

buildBWT: Makefile buildBWT.c $(BUILDBWT).c $(IO_LIST) $(MALLOC_COUNT_LIST) $(ITERATOR_LIST)
		$(CC) $(CFLAGS) $(BUILDBWT).c $(IO_OBJS) $(MALLOC_COUNT_OBJS) $(ITERATOR_OBJS) $(BWT_OBJS) $(LIBS) $(DIVSUFSORT_OBJS) -o $(BUILDBWT)


random:
mt19937ar.o: Makefile mt19937ar.c mt19937ar.h mt19937ar.o
	cd $(RANDOM_DIR) && $(CC) $(CFLAGS) $(LIBS) -c $(RANDOM_SRC)


# ---- CLEANING ----
clean:
	rm -f $(CALLBACKS_DIR)/*.o $(IO_DIR)/*.o $(ITERATOR_DIR)/*.o $(RANDOM_DIR)/*.o $(MALLOC_COUNT_DIR)/*.o $(PROGRAMS) 
