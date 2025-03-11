CC = gcc
CFLAGS = -lm -O2

# List all your source files here
SOURCES = A0_timestep_and_time.c \
			A1_combine_u_g_files.c \
			A3_get_RSL_sites.c

# Generate the names of the executables by removing the .c extension
EXECUTABLES = $(SOURCES:.c=)

# Default target
all: $(EXECUTABLES)

# Rule to compile each source file into an executable
%: %.c
	$(CC) $(CFLAGS) -o $@ $<

# Clean target
clean:
	rm -f $(EXECUTABLES)

.PHONY: all clean
