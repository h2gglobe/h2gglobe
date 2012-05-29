# Function to pull out a single value from a single tree from a single file and return it 
def getTreeEntry(fileName,treeName,branchName):

	try: ROOT.gROOT
	except NameError: import ROOT
	ROOT.gROOT.ProcessLine( \
	   "struct Entry{	\
    	    int r;	\
   	   };"
	)
	from ROOT import Entry
	newFile   = ROOT.TFile.Open(fileName)

	tree = newFile.Get(treeName)
	brnch= tree.GetBranch(branchName)
	valueStruct = Entry();
	brnch.SetAddress(ROOT.AddressOf(valueStruct,'r'))
	tree.GetEntry(0)

	newFile.Close()
	return int(valueStruct.r)

	
