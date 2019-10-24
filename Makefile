ROOTCONFIG	:= root-config
ObjSuf        = o
SrcSuf        = cpp
ExeSuf        =
DllSuf        = so
EVENTLIB      = $(EVENTO)
OutPutOpt     = -o
ROOTLIBS	:= $(shell $(ROOTCONFIG) --libs)
ROOTGLIBS	:= $(shell $(ROOTCONFIG) --glibs)
ROOTINC		:= $(shell $(ROOTCONFIG) --incdir)


INCDIRS = -I$(shell pwd) -I$(shell pwd) -I$(shell pwd)/BField -I${BOOST_INCLUDE}

# Linux
DBGFLAG       = -g
OPTFLAG       = -O2
CXX           = g++
CXXFLAGSDBG   = -g -Wall -fPIC -I$(ROOTINC) $(INCDIRS)
CXXFLAGSOPT   = -O2 -Wall -fPIC -fpermissive -I$(ROOTINC) $(INCDIRS)
LD            = g++
LDFLAGS       := $(shell $(ROOTCONFIG) --ldflags )
SOFLAGS       = -Wl,-soname,libEvent.so -shared 
LIBS          = $(ROOTLIBS) $(SYSLIBS) -lpython2.7 -lboost_system
GLIBS         = $(ROOTLIBS) $(ROOTGLIBS) $(SYSLIBS) 

CXXFLAGS      = $(CXXFLAGSOPT)

all: run runShared.so

default: run

PROGRAMS	= test

TESTO	= Reading.o test.o

FIELDSO = BField/BField.o BField/BList.o BField/BFieldGradient.o BField/BFieldROOT.o ParticleField.o BField/BDressing.o BField/cmsInterpnoiseGen.o BField/cmsB1PulseInterp.o BField/BFieldInterp.o

RUNO = $(FIELDSO) Neutron.o Boundary.o  Vector.o Reading.o Scattering.o Run.o RunParameters.o
RUNH = $(RUNO:.o=.h)
RUNDBG = $(RUNO:.o=.dbg)

tools: tools.o

test: $(TESTO)
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $@

run: bin/run

bin/run: $(RUNO) runDict.o run_main.o
	$(LD) -O2 $(LDFLAGS) $^ $(LIBS) -o $@

run_LinkDef.h: $(RUNH)
	python make_LinkDef.py -n run $(RUNH)

runShared.so: $(RUNO) runDict.o
	$(LD) --shared $(LDFLAGS) $^ $(LIBS) -o $@

run_dbg: $(RUNDBG)
	$(LD) -g $(LDFLAGS) $^ $(LIBS) -o bin/$@

clean:
	rm -f *.o *.dbg bin/run runDict.C runDict.h BField/*.o run_LinkDef.h


%.o:: %.C
	$(CXX) $(CXXFLAGS) -Wno-deprecated -c $< -o $@

%.dbg:: %.C
	$(CXX) $(CXXFLAGSDBG) -Wno-deprecated -c $< -o $@


runDict.C: $(RUNH) run_LinkDef.h
	@echo "Generating dictionary $@..."
	@rootcint -f $@ -c $^
