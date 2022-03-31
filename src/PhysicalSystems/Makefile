PHYS_SYS_OBJS = PhysicalSystemBase.o \
                Ising2D.o            \
                IsingND.o            \
                Heisenberg2D.o       \
                Heisenberg3D.o       \
                Elements.o           \
                CrystalBase.o        \
                CrystalStructure3D.o \
                Alloy3D.o            \
                HeisenbergHexagonal2D.o

.PHONY : default owl-qe clean 

default : owl

owl : libPhysicalSystems.a

owl-qe :
	cd QuantumEspresso && $(MAKE) owl-qe

clean :
	rm -rf *.o *.a *.dSYM
	@if test -f $(MASTER_DIR)/quantum_espresso/make.sys ; then	\
		echo ' Quantum Espresso exists. Cleaning up...' ;	\
		cd QuantumEspresso && $(MAKE) clean;			\
	 fi

%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE_PATH) -c -o $@ $<

libPhysicalSystems.a : $(PHYS_SYS_OBJS)
	$(AR) -r $@ $(PHYS_SYS_OBJS)
	@echo ' '
