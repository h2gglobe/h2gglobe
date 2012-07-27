CC        = g++
CCFLAGS   = -Wall
SOURCES   =
ROOTFLAGS = `root-config --cflags --glibs` -lTMVA


all:
	$(CC) $(CCFLAGS) -o TMVA_training $(ROOTFLAGS) TMVA_training.cc

TMVA_training:
	$(CC) $(CCFLAGS) -o TMVA_training $(ROOTFLAGS) TMVA_training.cc
clean:
	rm TMVA_training


#g++ -Wall -o TMVA_training `root-config --cflags --glibs`  -lTMVA TMVA_training.cc