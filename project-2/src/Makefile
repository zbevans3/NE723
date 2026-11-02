CXX      = g++
CXXFLAGS = -O3  -Wall -Wextra -MMD -MP

SRC    = main.cpp read_input.cpp build_mesh.cpp write_output.cpp \
         source_iteration.cpp outer_iteration.cpp
OBJ    = $(SRC:.cpp=.o)
DEP    = $(OBJ:.o=.d)
TARGET = project-2.exe

$(TARGET): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJ)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean debug
clean:
	@rm -f *.o *.d $(TARGET)

# Build with sanitizers to catch undefined behavior / memory errors
debug: CXXFLAGS = -O0 -g -Wall -Wextra -fsanitize=address,undefined -MMD -MP
debug: clean $(TARGET)

-include $(DEP)
