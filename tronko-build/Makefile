# the compiler: gcc for C program
CC = gcc

#compiler flags:
# -g adds debugging information to the executable file
# -Wall turns on most, but not all, compiler warnings
CFLAGS = -w -pg
DBGCFLAGS = -g -w
# -lm links the math library
#LIBS = -lm -lpthread -lz
LIBS = -lm -pthread -lz -lrt -std=gnu99
OPTIMIZATION = -march=native
#sources
SOURCES = tronko-build.c getclade.c readfasta.c readreference.c allocatetreememory.c math.c likelihood.c opt.c options.c printtree.c
HASHMAP = hashmap.c
OBJECTS = (SOURCES: .c = .o)
# the build target executable:
TARGET = tronko-build

all: $(TARGET)
$(TARGET): $(TARGET).c
	$(CC) $(OPTIMIZATION) -o $(TARGET) $(HASHMAP) $(SOURCES) $(LIBS)
debug: $(TARGET).c
	$(CC) $(DBGCFLAGS) -o $(TARGET) $(HASHMAP) $(SOURCES) $(LIBS)

clean:
	$(RM) $(TARGET)

