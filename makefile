CXX = icc
CXX_FLAGS = -O2 -std=c++14
EXEC = reaxdetect-dev
VPATH = reaxdetect reaxdetect/util reaxdetect/smiles 
OBJECTS = analyzer.o datawriter.o main.o reaxdetect.o reaxreader.o reaxreader_mol.o reaxreader_reac.o trajectory.o \
	smiles.o \
    path.o strutil.o

$(EXEC): $(OBJECTS)
	$(CXX) $(OBJECTS) -o $(EXEC)

REAXREADER_H = reaxreader.h trajectory.h typedef.h smiles.h simulation.h 
ANALYZER_H = analyzer.h $(REAXREADER_H)

analyzer.o : elements.h vecutil.h $(ANALYZER_H) 
datawriter.o : datawriter.h errors.h elements.h path.h strutil.h $(ANALYZER_H)
main.o : reaxdetect.h config.h $(ANALYZER_H)
reaxdetect.o : $(ANALYZER_H) datawriter.h errors.h reaxdetect.h elements.h path.h strutil.h
reaxreader.o : $(REAXREADER_H) algorithmutil.h path.h 
reaxreader_mol.o : $(REAXREADER_H) algorithmutil.h 
reaxreader_reac.o : $(REAXREADER_H) algorithmutil.h strutil.h 
trajectory.o : errors.h trajectory.h simulation.h algorithmutil.h strutil.h
smiles.o : smiles.h tablesort.h elementname.h
path.o : path.h
strutil.o : strutil.h

clean:
	rm -f $(EXEC) $(OBJECTS)
