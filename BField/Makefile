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


# Linux
DBGFLAG       = -g
OPTFLAG       = -O2
CXX           = g++
CXXFLAGSDBG   = -g -Wall -fPIC #-I$(ROOTINC)
CXXFLAGS      = -O2 -Wall -fPIC #-I$(ROOTINC)
LD            = g++
#LDFLAGS       := $(shell $(ROOTCONFIG) --ldflags)
SOFLAGS       = -Wl,-soname,libEvent.so -shared
#LIBS          = $(ROOTLIBS) $(SYSLIBS)
#GLIBS         = $(ROOTLIBS) $(ROOTGLIBS) $(SYSLIBS) 


default: run

PROGRAMS	= BField

BFIELD = BField.o BList.o

BDRESSING = $(BFIELD) BDressing.o

TESTO	= $(BFIELD)

test: $(BFIELD) test.o
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $@

dressing_py:
	echo -n '.L BField.C+\n.L BList.C+\n.L BDressing.C+\n' | root -l

FIELDSO = BField.o BList.o BDressing.o
FIELDSH = BField.h BList.h BDressing.h

fields.so: $(FIELDSO)
	$(LD) --shared $(LDFLAGS) $^ $(LIBS) -o $@

fieldsDict.C: $(FIELDSH) #fields_LinkDef.h
	rootcint -f $@ -c $^

clean:
	rm -f *.o *.dbg $(PROGRAMS)

%.o:: %.C
	$(CXX) $(CXXFLAGS) -Wno-deprecated -c $< -o $@

%.dbg:: %.C
	$(CXX) $(CXXFLAGSDBG) -Wno-deprecated -c $< -o $@


runDict.C: $(RUNH)
	@echo "Generating dictionary $@..."
	@rootcint -f $@ -c $^
