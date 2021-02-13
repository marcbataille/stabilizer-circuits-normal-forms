CC=gcc
CFLAGS= -Wall -Wconversion -Werror -Wextra -Wfatal-errors -Wpedantic -Wwrite-strings -O2
LDFLAGS=-lm
EXEC=stabnf

all: $(EXEC)

stabnf: main.o aux_ops.o cases.o print.o
	$(CC) -o $@ $^ $(LDFLAGS)

main.o: constants.h gate_struct.h print.h aux_ops.h cases.h

aux_ops.o: gate_struct.h

cases.o: gate_struct.h aux_ops.h

print.o: gate_struct.h constants.h aux_ops.h

%.o: %.c
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf *.o

mrproper: clean
	rm -rf $(EXEC)
