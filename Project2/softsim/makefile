CC=gcc
MAINDEPS = initializer.h running.h 
INITDEPS = globaldata.h initializer.h list.h timing.h
GLOBALDEPS = list.h
OBJS = main.o initializer.o running.o globaldata.o list.o timing.o
LIB = -lm
# -Wall

simprog : $(OBJS)
	$(CC) -g -o simprog $(OBJS) $(LIB)

main.o : main.c $(MAINDEPS) 
	$(CC) -g main.c -c 

initializer.o : initializer.c $(INITDEPS)
	$(CC) -g initializer.c -c 

running.o : running.c $(INITDEPS)
	$(CC) -g running.c -c 

globaldata.o : globaldata.c $(GLOBALDEPS)
	$(CC) -g globaldata.c -c 

list.o : list.c
	$(CC) -g list.c -c

time.o : timing.c $(GLOBALDEPS)
	$(CC) -g timing.c -c

.PHONY : clean
clean:
	rm -f ./*.o core *~

