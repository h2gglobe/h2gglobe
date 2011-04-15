HDR   = ./interface/
OBJ   = ./obj/
LIB   = ./lib/
BIN   = ./bin/
PRG   = ./test/
SRC   = ./src/
DIC   = $(SRC)

HdrSuf  =  .h
SrcSuf  =  .cc
ObjSuf  =  .o
PrgSuf  =  .cpp
BinSuf  =  .exe
LibSuf  =  .so

LIBNAME=libh2gglobeVertexAnalysis

HDRS  =  $(wildcard $(HDR)*$(HdrSuf))
SRCS  =  $(wildcard $(SRC)*$(SrcSuf))
## _OBJS  =  $(patsubst %$(SrcSuf), %$(ObjSuf), $(SRCS))
_OBJS  =  $(patsubst $(SRC)%$(SrcSuf), %$(ObjSuf), $(SRCS))
OBJS   =  $(addprefix $(OBJ),$(_OBJS))
PRGS  =  $(wildcard $(PRG)*$(PrgSuf))
_BINS    =  $(wildcard $(PRG)*$(PrgSuf))
__BINS   =  $(_BINS:$(PrgSuf)=$(BinSuf))
___BINS  =  $(notdir $(__BINS))
BINS	 =  $(addprefix $(BIN),${___BINS})

LINKDEF   =$(wildcard ${DIC}*LinkDef.h)
DICTHDRS  =$(LINKDEF)

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTGLIBS     = $(shell root-config --glibs)

CLHEPCFLAGS   =  $(shell clhep-config --include)
CLHEPLIBS     =  $(shell clhep-config --libs)

## to use a different TMVA wrt the one shipped with ROOT
TMVA_ROOT=$(HOME)/Analysis/tools/tmva/V04-01-00_root-5.22.00d-cms28
ifneq ($(TMVA_ROOT),)
TMVAFLAGS=-I$(TMVA_ROOT)
TMVALDD=-L$(TMVA_ROOT)/lib
endif

ifneq ($(DICTHDRS),)
OBJS+=$(OBJ)mydict.o
endif

ARCHL = -m32
ARCH  =  $(shell root-config --arch)
ifeq ($(ARCH),linuxx8664gcc)
ARCHL = -m64
endif

CXX  =  c++
CXXFLAGS  =  -Wall -O -fPIC -I$(HDR)  $(TMVAFLAGS) $(ROOTCFLAGS)
CPP  =  c++
CPPFLAGS  = -Wall $(ARCHL) -I$(HDR) $(TMVAFLAGS)  $(ROOTCFLAGS)

LD       =  c++
LDFLAGS  =  -rdynamic -shared $(ARCHL)
SONAME	 =  $(LIBNAME).so
SOFLAGS  =  -Wl,-soname,$(SONAME)

GLIBS   =  -lm -ldl -rdynamic $(TMVALDD) $(ROOTGLIBS) -lGenVector -lFoam -lMinuit -lTMVA -lMLP -lXMLIO  -lTreePlayer


#################################################
#if mac
ARCH  =  $(shell root-config --arch)
ifeq ($(ARCH),macosx64)
LibSuf  =  .dylib

CPPFLAGS  =  -Wall -W -Woverloaded-virtual -O2 $(ARCHL) -pipe -I$(HDR) $(ROOTCFLAGS)

CXXFLAGS  =  -Wall -W -Woverloaded-virtual -O2 $(ARCHL) -pipe -I$(HDR) $(ROOTCFLAGS)

LDFLAGS  =  -dynamiclib -shared -single_module -undefined dynamic_lookup $(ARCHL)
SONAME	 =  $(LIBNAME).dylib
SOFLAGS  =
endif
#################################################


.PHONY: all exe test clean


all: paths $(LIB)$(SONAME)


exe: paths $(BINS) 


test:
	@echo "HDRS = \"$(HDRS)\""
	@echo "DICTHDRS = \"$(DICTHDRS)\""
	@echo "SRCS = $(SRCS)"
	@echo "OBJS = $(OBJS)"
	@echo "PRGS = $(PRGS)"
	@echo "BINS = $(BINS)"

paths:
	mkdir -p $(LIB) $(OBJ)


$(OBJ)%$(ObjSuf): $(SRC)%$(SrcSuf) $(HDRS)
	$(CXX) -c $(CXXFLAGS) -o $@ $< 

$(OBJ)mydict.o: $(OBJ)mydict.cc
	$(CXX) -c $(CXXFLAGS) -I. -o $@ $<

$(OBJ)mydict.cc: $(DICTHDRS) $(HDRS)
	@echo "Generating dictionary for  ..."
	rootcint -f $(OBJ)mydict.cc -c -p ${CXXFLAGS} -I. $(HDRS) $(DICTHDRS)

$(LIB)$(SONAME): $(OBJS)
	@echo "Linking $(SONAME):"
	$(LD) $(LDFLAGS) $(OBJS) $(SOFLAGS) $(GLIBS) -o $(LIB)$(SONAME)

$(BIN)%$(BinSuf): $(PRG)%$(PrgSuf) $(HDRS) $(LIB)$(SONAME)
	$(CPP) $(CPPFLAGS) -L$(LIB) $(GLIBS) -lMyTest  -o $@ $<

clean:
	rm -f $(OBJ)*$(ObjSuf) $(LIB)*$(LibSuf) $(OBJ)mydict*
