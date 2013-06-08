#!/usr/bin/env python

fsm = open('table_sm.txt')
fgg = open('table_gg_grav.txt')
fqq = open('table_qq_grav.txt')

deets=[]

for f in [fsm,fgg,fqq]:
  info = {}
  for line in f.readlines():
    if line.startswith('-') or line.startswith('C'): continue
    name = line.split()[0]
    if name!='all': name=int(name.strip('cat'))
    info[name] = line.split()[1:]

  deets.append(info)

sm = deets[0]
gg = deets[1]
qq = deets[2]

arr = sm.keys()
arr.sort()
print arr
for key in arr:
  print '%6s'%key, '%5.2f'%(100.*float(sm[key][0])/float(sm['all'][0])), sm[key][5], sm[key][7],  '%5.2f'%(100.*float(gg[key][0])/float(gg['all'][0])), gg[key][5], gg[key][7], '%5.2f'%(100.*float(qq[key][0])/float(qq['all'][0])), qq[key][5], qq[key][7], '%5.1f'%float(sm[key][8]), sm[key][9], '%3.1f'%float(sm[key][10]), '%4.0f'%float(sm[key][12]) 
  
print '----------------------------------------------------------'
print ' for latex'
print '----------------------------------------------------------'

for key in arr:
  print '&', '%5.2f'%(100.*float(sm[key][0])/float(sm['all'][0])), '&', sm[key][5], '&', sm[key][7], '&', '%5.2f'%(100.*float(gg[key][0])/float(gg['all'][0])), '&', gg[key][5], '&', gg[key][7], '&', '%5.2f'%(100.*float(qq[key][0])/float(qq['all'][0])), '&', qq[key][5], '&', qq[key][7], '&', '%5.1f'%float(sm[key][8]), '&', '$\\pm$', '&', '%3.1f'%float(sm[key][10]), '&', '%4.0f'%float(sm[key][12]), '\\tabularnewline' 
    


