# Makefile enoise

CC = gcc
CFLAGS = -c
LDFLAGS = -lgsl -lgslcblas -lfftw3 -lrt -lm
SOURCES = enoise.c maker.c cmdparser.c stdata.c station.c tones.c
OBJECTS := $(SOURCES:.c=.o)
EXECUTABLE = enoise

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

clean:
	rm -rf *.o $(EXECUTABLE)
