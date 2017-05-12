CFLAGS:= -std=c99  -Wall
CC:=mpicc
LDLIBS:=  -lm 

INCL = -I ./

run: explicit_left_angle.o
	$(CC) $(CFLAGS)  -o $@ $< $(LDLIBS)
	@echo "--------------------------------"
	@echo "$@ is built successfully."


explicit_left_angle.o: explicit_left_angle.c
	$(CC) $(CFLAGS) $(LDLIBS) -c $<

clean:
	 rm -rf *.o *.d  run
