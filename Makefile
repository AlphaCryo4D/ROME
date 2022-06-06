#check the outside variable
ifeq (,${ROME_CC})
ROME_CC   := mpiicpc
endif
ifneq (,$(findstring mpi,${ROME_CC}))
MACROS    := -DUSEMPI
else
MACROS    :=
endif
#the complier and build option
CC	  := $(ROME_CC)
LD        := $(ROME_CC)

# build on KNL
ifeq (true,${ROME_KNL})
ifeq (true,${FLAT_MODE})
KNLFLAGS  := -xMIC-AVX512 -DUSEMCDRAM -I/opt/crtdc/memkind/latest/include/ /opt/crtdc/memkind/latest/lib/libmemkind.so
else
KNLFLAGS  := -xMIC-AVX512
endif
endif

# -O3 optimization may get some error for RELION(mask.cpp),NO PARALLEL MODE is used to compare with RELION
ifeq (true,${ROME_NO_PARALLEL})
FLAGS     := -mkl -parallel-source-info=2 -O2 -g $(TBB) -DNDEBUG -std=c++11 -static-intel $(KNLFLAGS)
else
FLAGS     := -fopenmp -mkl -parallel-source-info=2 -O2 -g -xHost $(TBB) -DNDEBUG -std=c++11 -static-intel $(KNLFLAGS)
endif

#debug for vtune and some vectorization report
ifeq (true,${ROME_DEBUG})
DEBUG     := -g -debug inline-debug-info -qopt-report-phase=vec -qopt-report=5
endif

#The following line disables offload.  By default we not offload, no flags needed
ifneq (true,${ROME_OFFLOAD})
OFFLOAD   := -qno-offload
endif
#This line allows offload debugging
#OFFLOAD   := -DSERIAL_OFFLOAD
TBB		  := 

SRC_DIR   := src
HEALPIX_DIR:= Healpix_2.15a
APPS_DIR  := apps#$(MODULES) #$(addprefix ./,$(MODULES))
RESMAP_DIR := src/resmap
BUILD_DIR := $(addprefix build/,$(SRC_DIR) $(APPS_DIR) $(HEALPIX_DIR) $(RESMAP_DIR))

# src files include map2d and gtm and other code.
SRC_CPP   := $(foreach sdir,$(SRC_DIR),$(wildcard $(sdir)/*.cpp))
HEALPIX_CC:= $(foreach sdir,$(HEALPIX_DIR),$(wildcard $(sdir)/*.cc))
MAP2D_CPP := $(foreach sdir,$(SRC_DIR),$(wildcard $(sdir)/map2d_*.cpp))
MAP3D_CPP := $(foreach sdir,$(SRC_DIR),$(wildcard $(sdir)/map3d_*.cpp))
RESMAP_CPP := $(foreach sdir,$(RESMAP_DIR), $(wildcard $(sdir)/*.cpp))
GTM_CPP   := $(foreach sdir,$(SRC_DIR),$(wildcard $(sdir)/gtm_*.cpp))
OTHER_CPP := $(filter-out $(MAP2D_CPP) $(MAP2D_OLD_CPP) $(MAP3D_CPP) $(GTM_CPP),$(SRC_CPP))

# rome_gtm,rome_map2d,rome_tool application 
APPS_CPP  := $(foreach sdir,$(APPS_DIR),$(wildcard $(sdir)/*.cpp))

SRC_OBJ   := $(patsubst %.cpp,build/%.o,$(SRC_CPP))
MAP2D_OBJ := $(patsubst %.cpp,build/%.o,$(MAP2D_CPP))
MAP3D_OBJ := $(patsubst %.cpp,build/%.o,$(MAP3D_CPP))
GTM_OBJ   := $(patsubst %.cpp,build/%.o,$(GTM_CPP))
RESMAP_OBJ := $(patsubst %.cpp,build/%.o,$(RESMAP_CPP))
OTHER_OBJ := $(patsubst %.cpp,build/%.o,$(OTHER_CPP))

HEALPIX_OBJ:=$(patsubst %.cc,build/%.o,$(HEALPIX_CC))

APPS_OBJ  := $(patsubst %.cpp,build/%.o,$(APPS_CPP))
APPS      := $(patsubst %.cpp,%,$(APPS_CPP))
APPS      := $(notdir $(APPS))

REOM_OBJ  := $(SRC_OBJ) $(HEALPIX_OBJ) 

ALL_APPS  := rome_map2d rome_sml rome_deep2d rome_tool rome_map3d rome_reconstruct rome_res

vpath %.cpp $(SRC_DIR) $(APPS_DIR) $(RESMAP_DIR)
vpath %.cc $(HEALPIX_DIR)

define make-goal
$1/%.o: %.cpp
	$(CC) $(FLAGS) $(MACROS) $(DEBUG) $(OFFLOAD) -c $$< -o $$@
$1/%.o: %.cc
	$(CC) $(FLAGS) $(MACROS) $(DEBUG) $(OFFLOAD) -c $$< -o $$@
endef

.PHONY: all checkdirs clean $(ALL_APPS)

all: $(ALL_APPS)

rome_res: checkdirs $(REOM_OBJ) $(RESMAP_OBJ) build/$(APPS_DIR)/rome_res.o
	$(LD) $(FLAGS) $(MACROS) $(OFFLOAD) build/$(APPS_DIR)/$@.o $(RESMAP_OBJ) $(HEALPIX_OBJ) -o bin/$@

rome_sml: checkdirs $(REOM_OBJ) build/$(APPS_DIR)/rome_sml.o
	$(LD) $(FLAGS) $(MACROS) $(OFFLOAD) build/$(APPS_DIR)/$@.o $(GTM_OBJ) $(OTHER_OBJ) $(HEALPIX_OBJ) -o bin/$@

rome_map2d: checkdirs $(REOM_OBJ) build/$(APPS_DIR)/rome_map2d.o
	$(LD) $(FLAGS) $(MACROS) $(OFFLOAD) build/$(APPS_DIR)/$@.o $(MAP2D_OBJ) $(OTHER_OBJ) $(HEALPIX_OBJ) -o bin/$@

rome_map3d: checkdirs $(REOM_OBJ) build/$(APPS_DIR)/rome_map3d.o
	$(LD) $(FLAGS) $(MACROS) $(OFFLOAD) build/$(APPS_DIR)/$@.o $(MAP3D_OBJ) $(OTHER_OBJ) $(HEALPIX_OBJ) -o bin/$@

rome_tool: checkdirs $(REOM_OBJ) build/$(APPS_DIR)/rome_tool.o
	$(LD) $(FLAGS) $(MACROS) $(OFFLOAD) build/$(APPS_DIR)/$@.o $(OTHER_OBJ) $(HEALPIX_OBJ) -o bin/$@

rome_deep2d: checkdirs $(REOM_OBJ) build/$(APPS_DIR)/rome_deep2d.o
	$(LD) $(FLAGS) $(MACROS) $(OFFLOAD) build/$(APPS_DIR)/$@.o $(MAP2D_OBJ) $(GTM_OBJ) $(OTHER_OBJ) $(HEALPIX_OBJ) -o bin/$@

rome_reconstruct: checkdirs $(REOM_OBJ) build/$(APPS_DIR)/rome_reconstruct.o
	$(LD) $(FLAGS) $(MACROS) $(OFFLOAD) build/$(APPS_DIR)/$@.o $(OTHER_OBJ) $(HEALPIX_OBJ) -o bin/$@

relion_refine: checkdirs $(RELION_NOMPI_OBJ) $(HEALPIX_OBJ) build/$(APPS_DIR)/relion_refine.o
	$(LD) $(FLAGS) $(MACROS) $(OFFLOAD) build/$(APPS_DIR)/$@.o $(RELION_NOMPI_OBJ) $(HEALPIX_OBJ) -o bin/$@

relion_refine_mpi: checkdirs $(RELION_OBJ) $(HEALPIX_OBJ) build/$(APPS_DIR)/relion_refine_mpi.o
	$(LD) $(FLAGS) $(MACROS) $(OFFLOAD) build/$(APPS_DIR)/$@.o $(RELION_OBJ) $(HEALPIX_OBJ) -o bin/$@

relion_reconstruct: checkdirs $(RELION_OBJ) $(HEALPIX_OBJ) build/$(APPS_DIR)/relion_reconstruct.o
	$(LD) $(FLAGS) $(MACROS) $(OFFLOAD) build/$(APPS_DIR)/$@.o $(RELION_OBJ) $(HEALPIX_OBJ) -o bin/$@

checkdirs: $(BUILD_DIR)

$(BUILD_DIR):
	@mkdir -p $@

clean:
	@rm -rf $(BUILD_DIR)
	@rm -rf ./build
	@rm ./bin/*
	
$(foreach bdir,$(BUILD_DIR),$(eval $(call make-goal,$(bdir))))
