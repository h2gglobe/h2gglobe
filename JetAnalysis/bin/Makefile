CC        = g++
CCFLAGS   = -Wall
SOURCES   =
ROOTFLAGS = `root-config --cflags --glibs` -lTMVA


all:
	$(CC) $(CCFLAGS) -o TMVA_training $(ROOTFLAGS) TMVA_training.cc
	$(CC) $(CCFLAGS) -o TMVA_Multiclass $(ROOTFLAGS) TMVA_Multiclass.cc
	$(CC) $(CCFLAGS) -o MakeJetIDPlots $(ROOTFLAGS) MakeJetIDPlots.cc JetTree.h

TMVA_training:
	$(CC) $(CCFLAGS) -o TMVA_training $(ROOTFLAGS) TMVA_training.cc

TMVA_Multiclass:
	$(CC) $(CCFLAGS) -o TMVA_Multiclass $(ROOTFLAGS) TMVA_Multiclass.cc

MakeJetIDPlots:
	$(CC) $(CCFLAGS) -o MakeJetIDPlots $(ROOTFLAGS) MakeJetIDPlots.cc JetTree.cc



clean:
	rm TMVA_training TMVA_Multiclass TMVA_application


#g++ -Wall -o TMVA_training `root-config --cflags --glibs`  -lTMVA TMVA_training.cc