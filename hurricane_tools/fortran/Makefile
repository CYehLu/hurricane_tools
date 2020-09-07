FILES = slp tk dbz vort interp

PREFIX = f90
TARGET = $(addprefix f90, $(FILES))


all : $(TARGET)

f90% : %.f90
	@echo ""
	@echo ""
	@echo "========================================================="
	@echo " Compiled '$<' as '$@' "
	@echo "========================================================="
	@echo ""
	f2py -c --opt='-O3' --f90flags='-fopenmp' -lgomp $< -m $@

clean :
	rm f90*.so
