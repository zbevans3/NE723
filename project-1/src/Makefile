CC = g++      # Compiler
CFLAGS=-O3 -Wall    # compiler options - turn on more error checking

SRC = main.cpp read_input.cpp build_mesh.cpp write_output.cpp source_iteration.cpp outer_iteration.cpp # Source files
HEADERS = data.h read_input.h build_mesh.h write_output.h source_iteration.h outer_iteration.h # Header files
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
build_mesh.o: data.h build_mesh.h
write_output.o: data.h write_output.h
source_iteration.o: data.h source_iteration.h
outer_iteration.o: data.h source_iteration.h outer_iteration.h