# instructions to produce a workspace from an opttree
. setup.sh
make
# You can take an example opttree from my public area.
# Otherwise a custom optree can be done running the code available in h2gglobe.
# If you split the files for different samples you have to hadd them and 
# finally you can run the merger.C macro (in UCSDPlotter) to produce the final tree
cp /afs/cern.ch/user/s/sani/public/opttree.root .
./plotter example_plotvariables.dat ./JimAnalysis.so opttree.root output.root

# The full configuration is in example_plotvariables.dat. 