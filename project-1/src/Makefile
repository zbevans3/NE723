CC = g++      # Compiler
CFLAGS=-O3 -Wall    # compiler options - turn on more error checking

SRC = main.cpp read_input.cpp  # Source files
HEADERS = data.h read_input.h  # Header files
TARGET = project-1.exe  # Output executable

OBJ=${SRC:.cpp=.o}

%.o: %.cpp
	$(CC) $(CFLAGS) -o $@ -c $<

$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ)

clean:
	@rm -f *.o $(TARGET)

# defne the header dependencies

main.o: data.h read_input.h 
read_input.o: data.h read_input.h
