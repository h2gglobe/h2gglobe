setrootdir () {
        export ROOTSYS=$1
        export PATH=$ROOTSYS/bin:$PATH
        export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
}

#----------------------------------------

export PATH="/afs/cern.ch/sw/lcg/external/Python/2.6.5/x86_64-slc5-gcc43-opt/bin:${PATH}"
export LD_LIBRARY_PATH="/afs/cern.ch/sw/lcg/external/Python/2.6.5/x86_64-slc5-gcc43-opt/lib:${LD_LIBRARY_PATH}"

# setup for libc/compiler
source /afs/cern.ch/sw/lcg/external/gcc/4.3.2/x86_64-slc5/setup.sh

# source /afs/cern.ch/sw/lcg/app/releases/ROOT/5.28.00g/x86_64-slc5-gcc43-opt/root/bin/thisroot.sh

setrootdir /afs/cern.ch/sw/lcg/app/releases/ROOT/5.28.00/x86_64-slc5-gcc43-opt/root
export PYTHONPATH=$ROOTSYS/lib
export ROOFITSYS=$ROOTSYS
export BOOST_DIR=/usr
