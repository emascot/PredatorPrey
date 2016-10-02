# Verbose
# make V=1
# Export
# make E=1
# Static load balancing
# make S=1

FC = mpif90
FFLAGS = -O3 -g -Jobj -Wall -cpp $(if $V,-DVERBOSE )$(if $E,-DEXPORT)$(if $S,-DSTATIC)
LIB =
FILES = parallel_tasks.f90 chase.f90
SRC = $(FILES:%=src/%)
OBJ = $(FILES:%.f90=obj/%.o)

all: chase_trajectory.out chase_time.out

chase_trajectory.out: src/chase_trajectory.f90 $(OBJ)
	$(FC) $(FFLAGS) $^ $(LIB) -o $@

chase_time.out: src/chase_time.f90 $(OBJ)
	$(FC) $(FFLAGS) $^ $(LIB) -o $@

$(OBJ): obj/%.o: src/%.f90 | obj
	$(FC) $(FFLAGS) -c $< -o $@

obj:
	mkdir -p obj

clean:
	rm -rf obj/* predprey.*
