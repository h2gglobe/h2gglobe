import os
import fnmatch
import sys

def printLinks(outFile,thispage,plotTypes,plotNames,wsName):
  outFile.write('<center>\n <p> \n <font size="5">Signal Interpolation Diagnostics</font>\n </p> \n')
  outFile.write('<a href=\"../../BMplots/grad/model.html\"><font size="4">View Background Model Diagnostics</font></a><br> \n')
  outFile.write('<script language="Javascript"> \n document.write("Results from interpolation run at: " + document.lastModified + ""); \n </SCRIPT> <br> \n ')
  outFile.write('Run on file: '+wsName+'<br>\n')
  outFile.write('Output file with interpolated histograms: '+wsName+'_interpolated.root <br>\n')
  for p, plot in enumerate(plotTypes):
    outFile.write('<a href=\"'+plot+'.html\">'+plotNames[p]+'</a><br> \n')
  outFile.write('<a href=\"../ada/'+thispage+'.html\">ada</a> \n')
  outFile.write('<a href=\"../grad/'+thispage+'.html\">grad</a> <br>\n')


wsName=sys.argv[1]
htmlFiles = []
paths=[]
plotTypes=[]

for root, dirs, files in os.walk('plots'):
  for dir in dirs:
    pathName=root+'/'+dir
    if pathName.count('/') < 2 and pathName not in paths: paths.append(pathName)
    elif dir not in plotTypes: plotTypes.append(dir)

plotNames=['Signal, background and data plots','Fractional systematic up and down templates','Up, down and interpolated signal plots','Background model in sidebands','Signal, background and data difference plots']

for path in paths:
  tempArr=path.split('/')
  bdtType=tempArr[1]

  for plot in plotTypes:
    allPlots=[]
    htmlFile = open(path+'/'+plot+'.html','w')
    printLinks(htmlFile,plot,plotTypes,plotNames,wsName)
    for root,dirs,files in os.walk(path+'/'+plot):
      for filename in fnmatch.filter(files,"*.png"):
        allPlots.append(filename)
    
    allPlots.sort()
    for p in allPlots:
      htmlFile.write('<a href='+plot+'/'+p+'><img height=\"400\" src=\"'+plot+'/'+p+'\"></a> \n')
