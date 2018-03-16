LIBS=`root-config --libs --cflags`
all:doFastGlauber.o FastGlauber.o Centrality.o CentralityMaker.o NegativeBinomial.o GlauberTree.o GlauberUtilities.o Nucleon.o
	@echo "compiling fastGlauber.exe ..." ;\
	g++ -o fastGlauber.exe doFastGlauber.o FastGlauber.o Centrality.o CentralityMaker.o NegativeBinomial.o 	GlauberTree.o GlauberUtilities.o Nucleon.o $(LIBS)

doFastGlauber.o:
	g++ -c doFastGlauber.cpp 

Nucleon.o:
	g++ -c FastGlauber/Nucleon.cxx $(LIBS)

FastGlauber.o:
	g++ -c FastGlauber/FastGlauber.cpp $(LIBS)

Centrality.o:
	g++ -c CentralityMaker/Centrality.cxx $(LIBS)

CentralityMaker.o:
	g++ -c CentralityMaker/CentralityMaker.cxx $(LIBS)

NegativeBinomial.o:
	g++ -c CentralityMaker/NegativeBinomial.cxx $(LIBS)

GlauberTree.o:
	g++ -c GlauberTree/GlauberTree.cxx $(LIBS)

GlauberUtilities.o:
	g++ -c GlauberUtilities/GlauberUtilities.cxx $(LIBS)


clean:
	@rm *.o ; \
	rm fastGlauber.exe
