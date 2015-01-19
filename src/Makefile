#
# Makefile for the extract_quantity tools
#
include ../../../models/Makefile.in
#
MODS=   $(DIR_COM)write_messages_mod.o \
	$(DIR_COM)physical_parameters_mod.o \
	$(DIR_COM)data_types_mod.o \
	$(DIR_COM)scene_creator_types_mod.o \
	$(DIR_COM)read_uff_mod.o \
	$(DIR_COM)interpol_sp.o \
	$(DIR_COM)read_scatt_info.o \
	$(DIR_COM)string_routines.o \
	$(DIR_COM)ref_ind_water.o \
	$(DIR_COM)read_master_lists_xml.o \
	$(DIR_COM)Ray_scatter_mod.o \
	$(DIR_COM)average_mod.o \
        $(DIR_COM)ncdf_utilities_mod.o \
	$(DIR_COM)fall_vel.o
#
MODS_nnc=   $(DIR_COM)write_messages_mod.o \
	$(DIR_COM)physical_parameters_mod.o \
	$(DIR_COM)data_types_mod.o \
	$(DIR_COM)scene_creator_types_mod.o \
	$(DIR_COM)read_uff_mod.o \
	$(DIR_COM)interpol_sp.o \
	$(DIR_COM)read_scatt_info.o \
	$(DIR_COM)string_routines.o \
	$(DIR_COM)ref_ind_water.o \
	$(DIR_COM)read_master_lists_xml.o \
	$(DIR_COM)Ray_scatter_mod.o \
	$(DIR_COM)average_mod.o 
#
OBJS=extract_quantity.o
OBJS2=extract_quantity_3d.o
OBJS3=L1_L2_rebin.o
OBJS4=L1_L2_compare.o
OBJS5=extract_quantity_hor.o
OBJS6=extract_env.o
OBJS7=uff_averager.o
OBJS8=uff_merger.o
#
$(MODS) :
	cd $(DIR_COM) && make $(MODS)

extract_quantity : $(MODS) $(OBJS)
	cd $(DIR_COM) && make $(MODS)
	$(FC) $(FFLAGS) $(MODS) $(OBJS) \
		-L$(XML_LIB) -L$(NC_LIB) \
		-lflib -lnetcdf -o ../bin/extract_quantity
#
extract_quantity.o : extract_quantity.f90
	cd $(DIR_COM) && make $(MODS)
	$(FC) $(FFLAGS) -I$(DIR_COM) -I$(XML_INCLUDE) -I$(NC_INCLUDE) -c extract_quantity.f90
#
extract_quantity_3d : $(MODS) $(OBJS2)
	cd $(DIR_COM) && make $(MODS)
	$(FC) $(FFLAGS) $(MODS) $(OBJS2) \
		-L$(XML_LIB) -L$(NC_LIB) \
		-lflib -lnetcdf -o ../bin/extract_quantity_3d
#
extract_quantity_3d.o : extract_quantity_3d.f90
	cd $(DIR_COM) && make $(MODS)
	$(FC) $(FFLAGS) -I$(DIR_COM) -I$(XML_INCLUDE) -I$(NC_INCLUDE) -c extract_quantity_3d.f90

#
extract_quantity_hor : $(MODS) $(OBJS5)
	cd $(DIR_COM) && make $(MODS)
	$(FC) $(FFLAGS) $(MODS) $(OBJS5) \
		-L$(XML_LIB) -L$(NC_LIB) \
		-lflib -lnetcdf -o ../bin/extract_quantity_hor
#
extract_quantity_hor.o : extract_quantity_hor.f90
	cd $(DIR_COM) && make $(MODS)
	$(FC) $(FFLAGS) -I$(DIR_COM) -I$(XML_INCLUDE) -I$(NC_INCLUDE) -c extract_quantity_hor.f90
#
L1_L2_rebin : $(MODS) $(OBJS3)
	cd $(DIR_COM) && make $(MODS)
	$(FC) $(FFLAGS) $(MODS) $(OBJS3) \
		-L$(XML_LIB) -L$(NC_LIB) \
		-lflib -lnetcdf -o ../bin/l1_l2_rebin
#
L1_L2_rebin.o : L1_L2_rebin.f90
	cd $(DIR_COM) && make $(MODS)
	$(FC) $(FFLAGS) -I$(DIR_COM) -I$(XML_INCLUDE) -I$(NC_INCLUDE) -c  L1_L2_rebin.f90
#
average_1_2_d.o : average_1_2_d.F90  
	cd $(DIR_COM) && make $(MODS)
	$(FC) $(FFLAGS) -I$(DIR_COM)  -I$(XML_INCLUDE)  -I$(NC_INCLUDE) -c  average_1_2_d.F90

L1_L2_compare : $(MODS) $(OBJS4)
	cd $(DIR_COM) && make $(MODS)
	$(FC) $(FFLAGS) $(MODS) $(OBJS4) \
		-L$(XML_LIB) -L$(NC_LIB) \
		-lflib -lnetcdf -o ../bin/l1_l2_compare
#
L1_L2_compare.o : L1_L2_compare.f90
	cd $(DIR_COM) && make $(MODS)
	$(FC) $(FFLAGS) -I$(DIR_COM) -I$(XML_INCLUDE) -I$(NC_INCLUDE) -c  L1_L2_compare.f90
#
extract_env : $(MODS) $(OBJS6)
	$(FC) $(FFLAGS) $(MODS) $(OBJS6) \
		-L$(XML_LIB) -L$(NC_LIB) \
		-lflib -lnetcdf -o ../bin/extract_env
#
extract_env.o : extract_env.f90
	cd $(DIR_COM) && make $(MODS)
	$(FC) $(FFLAGS) -I$(DIR_COM) -I$(XML_INCLUDE) -I$(NC_INCLUDE) -c extract_env.f90  
#
uff_averager : $(MODS_nnc) $(OBJS7)
	$(FC) $(FFLAGS) $(MODS_nnc) $(OBJS7) \
		-L$(XML_LIB) -lflib -o ../bin/uff_averager
#
uff_averager.o : uff_averager.f90
	cd $(DIR_COM) && make $(MODS_nnc)
	$(FC) $(FFLAGS) -I$(DIR_COM) -I$(XML_INCLUDE) -c uff_averager.f90
#
uff_merger : $(MODS) $(OBJS8)
	$(FC) $(FFLAGS) $(MODS_nnc) $(OBJS8) \
		-L$(XML_LIB) -lflib -o ../bin/uff_merger
#
uff_merger.o : uff_merger.f90
	cd $(DIR_COM) && make $(MODS_nnc)
	$(FC) $(FFLAGS) -I$(DIR_COM) -I$(XML_INCLUDE) -c uff_merger.f90
#
cleanup:
	mv $(DIR_COM)*.o   $(DIR_COM)../obj/
	mv $(DIR_COM)*.mod $(DIR_COM)../mod/
	mv *.o ../obj/
#
clean:
	rm -f $(DIR_COM)*.o   
	rm -f $(DIR_COM)../obj/*.o	
	rm -f $(DIR_COM)*.mod
	rm -f $(DIR_COM)../obj/*.o
	rm -f *.o  
	rm -f ../obj/*.o
