################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../DNA5_Basic_BWT.c \
../DNA5_tables.c \
../SLT.c \
../SLT_MAWs.c \
../dbwt.c \
../dbwt_queue.c \
../dbwt_utils.c \
../indexed_DNA5_seq.c \
../mt19937ar.c \
../naive_MAWs.c \
../sais.c \
../test_SLT_MAWs.c 

O_SRCS += \
../test_SLT_MAWs.o 

OBJS += \
./DNA5_Basic_BWT.o \
./DNA5_tables.o \
./SLT.o \
./SLT_MAWs.o \
./dbwt.o \
./dbwt_queue.o \
./dbwt_utils.o \
./indexed_DNA5_seq.o \
./mt19937ar.o \
./naive_MAWs.o \
./sais.o \
./test_SLT_MAWs.o 

C_DEPS += \
./DNA5_Basic_BWT.d \
./DNA5_tables.d \
./SLT.d \
./SLT_MAWs.d \
./dbwt.d \
./dbwt_queue.d \
./dbwt_utils.d \
./indexed_DNA5_seq.d \
./mt19937ar.d \
./naive_MAWs.d \
./sais.d \
./test_SLT_MAWs.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.c
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C Compiler'
	gcc -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


