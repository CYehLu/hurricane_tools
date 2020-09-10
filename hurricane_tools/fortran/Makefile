SRC = slp.f90 tk.f90 dbz.f90 vort.f90 interp.f90

NAMES = $(basename $(SRC))    # remove extension
TARGET = $(addprefix .f90, $(NAMES))


all : $(TARGET)

.f90% : %.f90
	@echo ""
	@echo ""
	@echo "========================================================="
	@echo " Compile '$<' as '$@' "
	@echo "========================================================="
	@echo ""
	rm -f -v $(@:.%=%)*.so      # remove old .so files
	f2py -c --opt='-O3' --f90flags='-fopenmp' -lgomp $< -m $@
	touch $@       # create tmp file to make sure that makefile can work well


.PHONY : clean
clean :
	rm -v f90*.so
	rm -v .f90*
