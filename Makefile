TARGET=condamage
HTS=../htslib
CFLAGS=-Wall -O2 -g -I$(HTS)
LDLIBS=-L$(HTS) -lm -Wl,-static -lhts -lz -Wl,-Bdynamic -pthread
CC=gcc

$(TARGET): $(TARGET).o
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

clean:
	rm -f $(TARGET) $(TARGET).o
