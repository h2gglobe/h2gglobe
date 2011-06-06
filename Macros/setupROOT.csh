


source /afs/cern.ch/sw/lcg/external/gcc/4.4.3/x86_64-slc5/setup.csh

setenv ROOTSYS /afs/cern.ch/sw/lcg/app/releases/ROOT/5.26.00b/x86_64-slc5-gcc44-opt/root

setenv PATH "$ROOTSYS/bin:$PATH"


if ($?LD_LIBRARY_PATH) then
   setenv LD_LIBRARY_PATH "$ROOTSYS/lib:$LD_LIBRARY_PATH"      # Linux, ELF HP-UX
else
   setenv LD_LIBRARY_PATH $ROOTSYS/lib
endif


if ($?PYTHONPATH) then
   setenv PYTHONPATH $ROOTSYS/lib:$PYTHONPATH
else
   setenv PYTHONPATH $ROOTSYS/lib
endif

#if ($?MANPATH) then
#   setenv MANPATH `dirname $ROOTSYS/man/man1`:$MANPATH
#else
#   setenv MANPATH `dirname $ROOTSYS/man/man1`
#endif



