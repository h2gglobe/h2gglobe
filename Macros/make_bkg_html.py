import os,sys,fnmatch

def printLinks(outFile,thispage,plotTypes,plotNames,wsName):
  outFile.write('<center> \n <p> \n <font size="5">Corrected Background Model Diagnostics</font> \n </p> \n')
  outFile.write('<a href=\"../../plots/grad/david.html\"><font size="4">View Signal Interpolation Diagnostics</font></a><br> \n')
  outFile.write('<script language="Javascript"> \n document.write("Results from correcting background model run at: " + document.lastModified + ""); \n </SCRIPT> <br> \n ')
  outFile.write('Rewritten workspace to file: '+wsName+'<br>\n')
  outFile.write('Output files: <br> \n')
  path = wsName[: wsName.find('CMS-HGG')]
  newWS = wsName[wsName.find('CMS-HGG') :]
  outFile.write(path+'bdtSidebandFits_ada_'+newWS+'<br>\n')
  outFile.write(path+'bdtSidebandFits_grad_'+newWS+'<br>\n')
  for p,plot in enumerate(plotTypes):
    outFile.write('<a href=\"'+plot+'.html\">'+plotNames[p]+'</a><br> \n')
  outFile.write('<a href=\"../ada/'+thispage+'.html\">ada</a>\n')
  outFile.write('<a href=\"../grad/'+thispage+'.html\">grad</a><br></center>\n')


wsName=sys.argv[1]

paths=[]
for root,dirs,files in os.walk('BMplots'):
  for dir in dirs:
    if root+'/'+dir not in paths: paths.append(root+'/'+dir)

oldMass="blank"

for path in paths:
  tempArr=path.split('/')
  bdtType=tempArr[1]
  plotTypes=['fCorr','fCovar','fit','uncorrErr','model']
  plotNames=['Fraction correlation matrices','Fractional covariance matrices','Bias Fits','Uncorrelated error matrices','Corrected background model']
  fCorrPlots=[]
  fCovarPlots=[]
  fitPlots=[]
  uncorrErrPlots=[]
  modelPlots=[]
  plots=[fCorrPlots,fCovarPlots,fitPlots,uncorrErrPlots,modelPlots]

  for root, dirs, files in os.walk(path):
    for p,pType in enumerate(plotTypes):
      for filename in fnmatch.filter(files,pType+"*.png"):
        plots[p].append(filename)
      
      plots[p].sort()
      htmlFile = open(path+'/'+plotTypes[p]+'.html','w')
      printLinks(htmlFile,plotTypes[p],plotTypes,plotNames,wsName);
      for plot in plots[p]:
        if plotTypes[p]=='fit':
          tArr=plot.split('_')
          t=tArr[1]
          mass=t[1:]
          if mass==oldMass:
            htmlFile.write('<a href='+plot+'><img width=\"300\" src=\"'+plot+'\"></a>\n')
          else:
            htmlFile.write('<br> <b> <font color=\"blue\"> Mass '+mass+':</font></b><br> <br> <a href='+plot+'><img width=\"300\" src=\"'+plot+'\"></a> \n')
          oldMass=mass
        else: 
          htmlFile.write('<a href='+plot+'><img height=\"400\" src=\"'+plot+'\"></a>\n')

  
