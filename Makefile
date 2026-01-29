# Makefile for MOND Physics Engine
# Supports multiple compilation targets and optimization levels

CC = gcc
CFLAGS = -Wall -Wextra -std=c11 -O3 -march=native -ffast-math
LDFLAGS = -lm
DEBUG_FLAGS = -g -O0 -DDEBUG
PROFILE_FLAGS = -pg

TARGET = mond_sim
SOURCE = mond_physics_engine.c

# Default target
all: $(TARGET)

# Standard optimized build
$(TARGET): $(SOURCE)
	$(CC) $(CFLAGS) -o $(TARGET) $(SOURCE) $(LDFLAGS)
	@echo "✓ Built optimized MOND simulator: $(TARGET)"

# Debug build with symbols
debug: $(SOURCE)
	$(CC) $(DEBUG_FLAGS) -o $(TARGET)_debug $(SOURCE) $(LDFLAGS)
	@echo "✓ Built debug version: $(TARGET)_debug"

# Profiling build
profile: $(SOURCE)
	$(CC) $(CFLAGS) $(PROFILE_FLAGS) -o $(TARGET)_profile $(SOURCE) $(LDFLAGS)
	@echo "✓ Built profiling version: $(TARGET)_profile"

# OpenMP parallel version (for future multi-threading)
parallel: $(SOURCE)
	$(CC) $(CFLAGS) -fopenmp -o $(TARGET)_parallel $(SOURCE) $(LDFLAGS)
	@echo "✓ Built parallel version: $(TARGET)_parallel"

# Run the simulation
run: $(TARGET)
	./$(TARGET)

# Run validation tests only
test: $(TARGET)
	@echo "Running MOND validation tests..."
	@./$(TARGET) < /dev/null

# Clean build artifacts
clean:
	rm -f $(TARGET) $(TARGET)_debug $(TARGET)_profile $(TARGET)_parallel
	rm -f *.o *.out *.dat
	rm -f gmon.out  # profiling output
	@echo "✓ Cleaned build artifacts"

# Clean everything including simulation output
cleanall: clean
	rm -f snapshot_*.dat rotation_curve*.dat
	@echo "✓ Cleaned all files"

# Install (copy to /usr/local/bin - requires sudo)
install: $(TARGET)
	sudo cp $(TARGET) /usr/local/bin/
	@echo "✓ Installed to /usr/local/bin/$(TARGET)"

# Show compiler info
info:
	@echo "Compiler: $(CC)"
	@echo "Flags: $(CFLAGS)"
	@echo "Libraries: $(LDFLAGS)"
	@$(CC) --version

# Benchmark (compile with timing and run)
benchmark: $(TARGET)
	@echo "Running benchmark..."
	time ./$(TARGET) < /dev/null

# Help message
help:
	@echo "MOND Physics Engine - Build System"
	@echo ""
	@echo "Available targets:"
	@echo "  make           - Build optimized version"
	@echo "  make debug     - Build debug version with symbols"
	@echo "  make profile   - Build with profiling support"
	@echo "  make parallel  - Build with OpenMP support"
	@echo "  make run       - Build and run simulation"
	@echo "  make test      - Run validation tests"
	@echo "  make benchmark - Run performance benchmark"
	@echo "  make clean     - Remove build artifacts"
	@echo "  make cleanall  - Remove all generated files"
	@echo "  make install   - Install to system (requires sudo)"
	@echo "  make help      - Show this message"

.PHONY: all debug profile parallel run test clean cleanall install info benchmark help
