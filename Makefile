DIR_MAIN       = ./
DIR_SRC        = $(DIR_MAIN)build/src
DIR_H          = $(DIR_MAIN)build/include
DIR_BUILD      = $(DIR_MAIN)build
DIR_OBJ        = $(DIR_BUILD)/obj

DEBUG = 
OPTIMIZATION = -O2
FLOWTRACE =
OPTIONS = -std=c++17
LINK_OPTIONS = 
CFLAGS = $(DEBUG) $(OPTIMIZATION) $(FLOWTRACE) $(OPTIONS)
COMPILER = g++
LIBS = -fopenmp -lgsl -lgslcblas -lm -lcuba
INCLUDES = -I $(DIR_H) 

CPP := $(shell find $(DIR_SRC) -name '*.cpp' -not -name '.*')
OBJ  = $(CPP:$(DIR_SRC)%.cpp=$(DIR_OBJ)%.o)

EXE =\
quenching

$(EXE): $(OBJ)
	echo "Linking:   $@ ($(COMPILER))"
	$(COMPILER) $(LINK_OPTIONS) -o $@ $^ $(LIBS) $(INCLUDES) 

$(DIR_OBJ)%.o: $(DIR_SRC)%.cpp
	@[ -d $(DIR_OBJ) ] || mkdir -p $(DIR_OBJ)
	@echo "Compiling: $< ($(COMPILER))"
	$(COMPILER) $(CFLAGS) $(INCLUDES) -c -o $@ $<

$(DIR_OBJ)%.o: $(DIR_SRC)%.cu
	@[ -d $(DIR_OBJ) ] || mkdir -p $(DIR_OBJ) 
	@echo "Compiling: $< ($(COMPILER))"
	$(COMPILER) $(CFLAGS) $(INCLUDES) -c -o $@ $<

clean:
	@echo "Object files and executable deleted"
	if [ -d "$(DIR_OBJ)" ]; then rm -rf $(EXE) $(DIR_OBJ)/*; rmdir $(DIR_OBJ); fi
	if [ -d "$(DIR_MAIN)/output" ]; then \
		if [ -z "$(shell find $(DIR_MAIN)/output -mindepth 1)" ]; then \
			rm -rf $(DIR_MAIN)/output; \
		else \
			echo "$(DIR_MAIN)/output is not empty. Do you want to delete it? (y/n)"; \
			read ans; \
			if [ "$$ans" = "y" ]; then rm -rf $(DIR_MAIN)/output; fi; \
		fi; \
	fi

.SILENT:
