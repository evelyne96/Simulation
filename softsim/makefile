CC=gcc
MAINDEPS = initializer.h running.h 
INITDEPS = globaldata.h initializer.h list.h
GLOBALDEPS = list.h
OBJS = main.o initializer.o running.o globaldata.o list.o
LIB = -lm
# -Wall

simprog : $(OBJS)
	$(CC) -o simprog $(OBJS) $(LIB)

main.o : main.c $(MAINDEPS) 
	$(CC) main.c -c 

initializer.o : initializer.c $(INITDEPS)
	$(CC) initializer.c -c 

running.o : running.c $(INITDEPS)
	$(CC) running.c -c 

globaldata.o : globaldata.c $(GLOBALDEPS)
	$(CC) globaldata.c -c 

list.o : list.c
	$(CC) list.c -c

.PHONY : clean
clean:
	rm -f ./*.o core *~

