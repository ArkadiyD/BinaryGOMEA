CPP := g++

# Folders
SRC_DIR := src
INC_DIR := include
BUILD_DIR := build
TARGET_DIR := bin

TARGET := GOMEA

SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(BUILD_DIR)/%.o,$(SRC_FILES))
$(info $$found .cpp files: [${SRC_FILES}])
$(info $$expected .o files: [${OBJ_FILES}])

CPPFLAGS := -std=c++17 -O3  
INC := -I$(INC_DIR)
LIBFLAGS := -lstdc++fs

$(TARGET): $(OBJ_FILES)
	$(CPP) -o $@ $^ $(LIBFLAGS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CPP) $(CPPFLAGS) $(INC) -c -o $@ $< $(LIBFLAGS) 

clean:
	@echo "Cleaning..." 
	$(RM) -r $(TARGET); $(RM) -r $(BUILD_DIR)/*.o 

.PHONY: clean