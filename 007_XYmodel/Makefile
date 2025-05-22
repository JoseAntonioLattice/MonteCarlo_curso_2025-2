
program=1d_xy_model

source= precision pbc statistics id_cluster constants random parameters spin_field measurements update_algorithms dynamics main
objects= $(patsubst %, bin/%.o,$(source))

bin/$(program): $(objects)
	gfortran $^ -o $@

bin/%.o: src/%.f90
	gfortran -c -J bin $< -o $@

bin/main.o: app/main.f90
	gfortran -c -J bin $< -o $@

run :
	echo "input/input_parameters.nml" | bin/$(program)

clean :
	rm bin/*
