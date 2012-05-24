import os
import sys

filename=sys.argv[1]
outputdir=sys.argv[2]

os.system("cp mva-plots-ada/* BMplots/ada")
os.system("cp mva-plots-grad/* BMplots/grad")
os.system("cp FitPlots/* BMplots/ada")
os.system("cp FitPlots/* BMplots/grad")

os.system("python make_html.py "+filename)
os.system("python make_bkg_html.py "+filename)

if not os.path.isdir(outputdir):
  os.makedirs(outputdir)

os.system("cp -r plots "+outputdir)
os.system("cp -r BMplots "+outputdir)

