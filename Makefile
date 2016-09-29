# Verbose
# make V=1
# Export
# make E=1

FC = mpif90
FFLAGS = -O3 -g -Jobj -Wall -cpp $(if $V,-DVERBOSE )$(if $E,-DEXPORT)
LIB =
FILES = parallel_tasks.f90 main.f90
SRC = $(FILES:%=src/%)
OBJ = $(FILES:%.f90=obj/%.o)

predprey.out: $(OBJ)
	$(FC) $(FFLAGS) $^ $(LIB) -o $@

$(OBJ): obj/%.o: src/%.f90 | obj
	$(FC) $(FFLAGS) -c $< -o $@

obj:
	mkdir -p obj

clean:
	rm -rf obj/* predprey.*
