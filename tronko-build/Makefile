CC = gcc

# ==== SOURCE FILES ====
SOURCES = tronko-build.c \
          getclade.c \
          readfasta.c \
          readreference.c \
          allocatetreememory.c \
          mathtools.c \
          likelihood.c \
          opt.c \
          options.c \
          printtree.c \
          hashmap/hashmap.c

TARGET  = tronko-build

# ==== BUILD DIRECTORIES ====
BUILD_DIR       = build
RELEASE_DIR     = $(BUILD_DIR)/release
DEBUG_DIR       = $(BUILD_DIR)/debug

# ==== FLAGS ====
# Warnings/standard shared by both configurations
WARN_FLAGS   = -Wall -Wextra -Wshadow -Wno-unused-result -Wno-format-truncation -Wno-format-overflow
STD_FLAGS    = -std=gnu99

# Optimization / debug flags
OPT_FLAGS    = -O3
DBG_FLAGS    = -g -O0

# Preprocessor flags: -MMD -MP auto-generate .d dependency files for headers
CPPFLAGS     = -MMD -MP

# Libraries
LDLIBS       = -lm -pthread -lz -lgsl -lgslcblas

CFLAGS_RELEASE = $(WARN_FLAGS) $(STD_FLAGS) $(OPT_FLAGS)
CFLAGS_DEBUG   = $(WARN_FLAGS) $(STD_FLAGS) $(DBG_FLAGS)

# ==== OBJECT FILES ====
OBJECTS_RELEASE = $(addprefix $(RELEASE_DIR)/, $(SOURCES:.c=.o))
OBJECTS_DEBUG   = $(addprefix $(DEBUG_DIR)/,   $(SOURCES:.c=.o))

# Dependency files generated alongside objects
DEPS_RELEASE = $(OBJECTS_RELEASE:.o=.d)
DEPS_DEBUG   = $(OBJECTS_DEBUG:.o=.d)

.PHONY: all debug clean

# ==== TARGETS ====
all: $(TARGET)

$(TARGET): $(OBJECTS_RELEASE)
	$(CC) $(CFLAGS_RELEASE) -o $@ $^ $(LDLIBS)

debug: $(TARGET)_debug

$(TARGET)_debug: $(OBJECTS_DEBUG)
	$(CC) $(CFLAGS_DEBUG) -o $@ $^ $(LDLIBS)

$(RELEASE_DIR)/%.o: %.c | $(RELEASE_DIR)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS_RELEASE) $(CPPFLAGS) -c $< -o $@

$(DEBUG_DIR)/%.o: %.c | $(DEBUG_DIR)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS_DEBUG) $(CPPFLAGS) -c $< -o $@

$(RELEASE_DIR) $(DEBUG_DIR):
	mkdir -p $@

# Pull in auto-generated header dependencies (silently ignored if missing)
-include $(DEPS_RELEASE)
-include $(DEPS_DEBUG)

clean:
	$(RM) -r $(BUILD_DIR) $(TARGET) $(TARGET)_debug
