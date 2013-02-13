import sys

nCats=4
errFile = open(sys.argv[1])

totEffs=[]
totEffErrs=[]

for line in errFile.readlines():
  eff=[]
  effErr=[]
  if line.startswith('#'): continue
  if len(line.split())==0: continue
  cat = line.split()[0]
  for i,el in enumerate(line.split()):
    if i==0: continue
    if (i%2==0): effErr.append(float(el))
    else: eff.append(float(el))

  #print cat, 'eff:', eff
  #print cat, 'err:', effErr
  
  if len(eff)!=len(effErr): sys.exit('ERROR: arrays not equal length')
  
  toteff = 1.;
  toteffErrSq = 0.;
  for i, e in enumerate(eff):
    toteff*=e
    toteffErrSq+=((effErr[i]/e)**2)

  toterr = toteff*(toteffErrSq**0.5)
  
  totEffs.append(toteff)
  totEffErrs.append(toterr)

  #print 'Total eff %s: %1.3f +/- %1.3f'%(cat,toteff,toterr)

print '\nFor mergeGlobeSystematics.C :'
print 'Double_t ratioTP_[nphocats]            = {%1.3f,%1.3f,%1.3f,%1.3f};'%(totEffs[0],totEffs[1],totEffs[2],totEffs[3])
print 'Double_t ratioTP_low_err_[nphocats]    = {%1.3f,%1.3f,%1.3f,%1.3f};'%(totEffErrs[0],totEffErrs[1],totEffErrs[2],totEffErrs[3])
print 'Double_t ratioTP_high_err_[nphocats]   = {%1.3f,%1.3f,%1.3f,%1.3f};'%(totEffErrs[0],totEffErrs[1],totEffErrs[2],totEffErrs[3])
