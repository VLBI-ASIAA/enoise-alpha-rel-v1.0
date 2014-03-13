# Makefile genran
# Author: Zheng Meyer-Zhao
# Date: 2013-06-27

CC = gcc
CFLAGS = -c
LDFLAGS = -lgsl -lgslcblas
SOURCES = genran.c
OBJECTS := $(SOURCES:.c=.o)
EXECUTABLE = genran

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

%.o: %.c
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@

clean:
	rm -rf *.o $(EXECUTABLE)
#genran: genran.o
#	cc genran.o -o genran -lgsl -lgslcblas
#genran.o:
#	cc -c genran.c
