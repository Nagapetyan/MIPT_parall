all: task1 task2 task3

CFLAGS:= -std=c99  -Wall
CC:=mpicc
LDLIBS:= -lgmp -lm 

INCL = -I ./

task3: integral.o
	$(CC) $(CFLAGS)  -o  $@ $< 
	@echo "--------------------------------"
	@echo "$@ is built successfully."

task2: task2.o
	$(CC) $(CFLAGS)  -o  $@ $< $(LDLIBS)
	@echo "--------------------------------"
	@echo "$@ is built successfully."


task1: task1.o
	$(CC) $(CFLAGS)  -o $@ $< $(LDLIBS)
	@echo "--------------------------------"
	@echo "$@ is built successfully."

integral.o: integral.c
	$(CC) $(CFLAGS)  -c $<

task2.o: task2.c
	$(CC) $(CFLAGS) $(LDLIBS) -c $<


task1.o: task1.c
	$(CC) $(CFLAGS) $(LDLIBS) -c $<

clean:
	 rm -rf *.o *.d  

.PHONY: clean