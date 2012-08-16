CC        = g++
CCFLAGS   = -Wall
SOURCES   =
ROOTFLAGS = `root-config --cflags --glibs` -lTMVA


all:
	$(CC) $(CCFLAGS) -o TMVA_training $(ROOTFLAGS) TMVA_training.cc
	$(CC) $(CCFLAGS) -o TMVA_Multiclass $(ROOTFLAGS) TMVA_Multiclass.cc

TMVA_training:
	$(CC) $(CCFLAGS) -o TMVA_training $(ROOTFLAGS) TMVA_training.cc

TMVA_Multiclass:
	$(CC) $(CCFLAGS) -o TMVA_Multiclass $(ROOTFLAGS) TMVA_Multiclass.cc

clean:
	rm TMVA_training TMVA_Multiclass


#g++ -Wall -o TMVA_training `root-config --cflags --glibs`  -lTMVA TMVA_training.cc