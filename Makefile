#
#  This file is part of synth.
#
#  synth is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  synth is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with synth.  If not, see <http://www.gnu.org/licenses/>.
#


OB_INC=../openbabel-2.3.2/include/
OB_LIB_DIR=.
OB_LIB=openbabel
#GSL_DIR=/root/gsl-1.16/
GSL_LIB=gsl
CBLAS_LIB=gslcblas
ZLIB=z

IDIR =./
CC=g++
OPT= -O2 -g # -p 
CFLAGS= $(OPT) -I$(IDIR) -I$(OB_INC) -L$(OB_LIB_DIR) -l$(ZLIB) -l$(GSL_LIB) -l$(CBLAS_LIB) -l$(OB_LIB) -lpthread
ODIR=./obj

DEPS = EdgeAggregator.h \
	EdgeAnnotation.h \
	Instantiator.h \
	Linker.h \
	Molecule.h \
	Rigid.h \
	Utilities.h \
	Atom.h \
	ConnectableAtom.h \
	LinkerConnectableAtom.h \
	RigidConnectableAtom.h \
	Bond.h \
	AtomT.h \
	IdFactory.h \
	OBWriter.h \
	Options.h \
	Validator.h \
	obgen.h \
	Constants.h \
	Thread_Pool.h \
	LevelHashMap.h \
	LikeMoleculesContainer.h \
	MinimalMolecule.h \
	SmiMinimalMolecule.h \
	MoleculeHashHypergraph.h \
	FixedSortedList.h \
	EdgeDatabase.h \
	FragmentEdgeMap.h \
	SimpleFragmentGraph.h \
	zpipe.h \
	TimedHashMap.h \
	TimedLikeValueContainer.h \
	bloom_filter.hpp

_OBGEN_DEPS = obgen.h 

_WRITER_DEPS = compliantwriter.h 

_TOSMI_DEPS = convertToSMI.h 

_OBJ = Atom.o \
	Instantiator.o \
	Linker.o\
	Main.o \
	Molecule.o \
	Rigid.o \
	Utilities.o \
	Atom.o \
	ConnectableAtom.o \
	LinkerConnectableAtom.o \
	RigidConnectableAtom.o \
	Bond.o \
	AtomT.o \
	OBWriter.o \
	IdFactory.o \
	Options.o \
	Validator.o \
	obgen.o \
	Constants.o \
	FragmentEdgeMap.o \
	zpipe.o \
	TimedHashMap.o \
	TimedLikeValueContainer.o


OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

synth: $(OBJ)
	$(CC) $^ $(CFLAGS) -o $@



_OBGEN_OBJ = synthobgen.o

OBGEN_OBJ = $(patsubst %,$(ODIR)/%,$(_OBGEN_OBJ))

$(ODIR)/%.o: %.cpp $(_OBGEN_DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

synthobgen: $(OBGEN_OBJ)
	$(CC) $^ $(CFLAGS) -o $@



_WRITER_OBJ = compliantwriter.o

WRITER_OBJ = $(patsubst %,$(ODIR)/%,$(_WRITER_OBJ))

$(ODIR)/%.o: %.cpp $(_WRITER_DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

compliantwriter: $(WRITER_OBJ)
	$(CC) $^ $(CFLAGS) -o $@




_TOSMI_OBJ = convertToSMI.o

TOSMI_OBJ = $(patsubst %,$(ODIR)/%,$(_TOSMI_OBJ))

$(ODIR)/%.o: %.cpp $(_TOSMI_DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

convertToSMI: $(TOSMI_OBJ)
	$(CC) $^ $(CFLAGS) -o $@



.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core synth.exe synth.exe.stackdump $(INCDIR)/*~
