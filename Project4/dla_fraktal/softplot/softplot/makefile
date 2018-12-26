CC=g++
CFLAGS=-I.
MAINDEPS = globaldata.h drawing.h 
DRAWDEPS = globaldata.h
OBJS = main.o drawing.o globaldata.o
LIB = -lm -lglut -lGLU -lGL
# -Wall

plotprog : $(OBJS)
	$(CC) -o plotprog $(OBJS) $(LIB)

main.o : main.cpp $(MAINDEPS) 
	$(CC) main.cpp -c 

drawing.o : drawing.cpp $(INITDRAWDEPSDEPS)
	$(CC) drawing.cpp -c 

globaldata.o : globaldata.cpp 
	$(CC) globaldata.cpp -c 

.PHONY : clean
clean:
	rm -f ./*.o core *~

