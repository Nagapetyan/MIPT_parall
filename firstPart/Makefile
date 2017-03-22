all: task1 task2

CFLAGS:= -std=c99  -Wall
CC:=mpicc
LDLIBS:= -lgmp -lm 

INCL = -I ./


task2: task2.o
	$(CC) $(CFLAGS)  -o  $@ $< $(LDLIBS)
	@echo "--------------------------------"
	@echo "$@ is built successfully."


task1: task1.o
	$(CC) $(CFLAGS)  -o $@ $< $(LDLIBS)
	@echo "--------------------------------"
	@echo "$@ is built successfully."


task3s.o: task2.c
	$(CC) $(CFLAGS) $(LDLIBS) -c $<


task1.o: task1.c
	$(CC) $(CFLAGS) $(LDLIBS) -c $<

clean:
	 rm -rf *.o *.d  

.PHONY: clean