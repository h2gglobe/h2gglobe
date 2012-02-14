#Building Signal Model

#need root and RooFit
#In CMSSW area, just do 
cmsenv

#On lxplus do 
export BOOST_DIR=/usr

#Make the area:
make


# Edit your copy of example-config.xml for your workspace



./fitSignal my-config.xml
