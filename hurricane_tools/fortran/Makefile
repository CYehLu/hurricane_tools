SRC1 = slp.f90 tk.f90 dbz.f90 vort.f90 interpz.f90
SRC2 = module_xy2rt.f90 module_rt2xy.f90

TARGET1 = $(SRC1:%.f90=.f90%)
TARGET2 = $(SRC2:module_%.f90=.f90%)


all : $(TARGET1) $(TARGET2)

.f90xy2rt : module_xy2rt.f90
	@echo ""
	@echo ""
	@echo "========================================================="
	@echo " Compile '$<' as '$@' "
	@echo "========================================================="
	@echo ""
	@echo ""
	rm -f -v $(@:.%=%)*.so     # remove old .so files
	f2py -c --opt='-O3' $< -m $(@:.%=%)
	touch $@       # create temporary files to make sure that makefile can work well

.f90rt2xy : module_rt2xy.f90
	@echo ""
	@echo ""
	@echo "========================================================="
	@echo " Compile '$<' as '$@' "
	@echo "========================================================="
	@echo ""
	@echo ""
	rm -f -v $(@:.%=%)*.so    
	f2py -c --opt='-O3' $< -m $(@:.%=%)
	touch $@     

.f90% : %.f90
	@echo ""
	@echo ""
	@echo "========================================================="
	@echo " Compile '$<' as '$@' "
	@echo "========================================================="
	@echo ""
	rm -f -v $(@:.%=%)*.so     
	f2py -c --opt='-O3' --f90flags='-fopenmp' -lgomp $< -m $(@:.%=%)
	touch $@     


.PHONY : clean
clean :
	rm -v f90*.so
	rm -v .f90*
