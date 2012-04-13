
source /afs/cern.ch/sw/lcg/external/gcc/4.4.3/x86_64-slc5/setup.sh
export ROOTSYS=/afs/cern.ch/sw/lcg/app/releases/ROOT/5.26.00b/x86_64-slc5-gcc44-opt/root
export PATH="$ROOTSYS/bin:$PATH"


if [ -n "${LD_LIBRARY_PATH:+x}" ]
then
   export LD_LIBRARY_PATH="$ROOTSYS/lib:$LD_LIBRARY_PATH"      # Linux, ELF HP-UX
else
   export LD_LIBRARY_PATH=$ROOTSYS/lib
fi


if [ -n "${PYTHONPATH:+x}" ] 
then
   export PYTHONPATH=$ROOTSYS/lib:$PYTHONPATH
else
   export PYTHONPATH=$ROOTSYS/lib
fi

#if ($?MANPATH) then
#   setenv MANPATH `dirname $ROOTSYS/man/man1`:$MANPATH
#else
#   setenv MANPATH `dirname $ROOTSYS/man/man1`
#endif



