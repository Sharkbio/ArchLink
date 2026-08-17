CC ?=	cc
CFLAG= -O3
SRCS=	util_lib.c hmm_lib.c run_hmm.c
OBJ=	util_lib.o hmm_lib.o run_hmm.o 

all:	fgs

fgs:	$(OBJ)
	$(CC)  $(CFLAG) -o FragGeneScan util_lib.o hmm_lib.o run_hmm.o  -lm -lpthread

%.o:	%.c
	$(CC) $(CFLAG) -c $< -o $@

clean:
	rm -rf *.o FragGeneScan* *~
