CFLAGS:= -std=c99  -Wall
CC:=mpicc
LDLIBS:=  -lm 

INCL = -I ./

run: merge_sort.o
	$(CC) $(CFLAGS)  -o $@ $< $(LDLIBS)
	@echo "--------------------------------"
	@echo "$@ is built successfully."


merge_sort.o: merge_sort.c
	$(CC) $(CFLAGS) $(LDLIBS) -c $<

clean:
	 rm -rf *.o *.d  run

.PHONY: clean