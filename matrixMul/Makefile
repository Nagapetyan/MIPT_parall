CFLAGS:= -std=c99  -Wall
CC:=mpicc
LDLIBS:= -lgmp -lm 

INCL = -I ./

run: mul.o
	$(CC) $(CFLAGS)  -o $@ $< $(LDLIBS)
	@echo "--------------------------------"
	@echo "$@ is built successfully."


mul.o: mul.c
	$(CC) $(CFLAGS) $(LDLIBS) -c $<

clean:
	 rm -rf *.o *.d  

.PHONY: clean