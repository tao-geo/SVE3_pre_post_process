CC = gcc
CFLAGS = -lm

# List all your source files here
SOURCES = A1_epoch_to_timestep.c \
		A2_get_RSL_c.c \
		combine_u_g_files.c \
		A3_refine_ocnfunc.c

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
