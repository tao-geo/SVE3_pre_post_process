CC = gcc
TARGET = precalc_oceanfunc

all: $(TARGET).x

$(TARGET).x: $(TARGET).c
	$(CC) -o $@ $< -lm

clean:
	rm -f $(TARGET)
