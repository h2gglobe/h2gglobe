import os
import sys

if len(sys.argv)!=4:
	sys.exit('usage mk_spin_cat_cards.py <card> <kinCats> <spinCats>')

def checkMulti(cardname,inc_cats):
	os.system('mv %s temp.txt'%cardname)
	newf = open(cardname,'w')
	oldf = open('temp.txt')
	for line in oldf.readlines():
		if 'discrete' in line:
			for cat in inc_cats:
				if 'pdfindex_%s'%(cat.split('cat')[1]) in line: newf.write(line)
				else: continue
		else:
			newf.write(line)
	newf.close()
	oldf.close()
	os.system('rm -f temp.txt')

card = open(sys.argv[1])
kinCats = int(sys.argv[2])
spinCats = int(sys.argv[3])

cats = set()
lines = card.readlines()
for line in lines:
	if line.startswith('bin'):
		for el in line.split()[1:]:
			cats.add(el)
		break

print cats

for sCat in range(spinCats):
	exc_string=''
	inc_cats=[]
	for kCat in range(kinCats):
		cat = kCat*spinCats+sCat
		for gcat in cats:
			if 'cat%d_'%cat in gcat:
				inc_cats.append(gcat)

	print inc_cats
	for cat in cats:
		if cat not in inc_cats:
			exc_string += 'ch1_%s|'%cat
	exc_string = exc_string[:-1]
	newcardname = card.name.replace('.txt','_spinCat%d.txt'%sCat)
	os.system('combineCards.py %s --xc=\"%s\" > %s'%(card.name,exc_string,newcardname))
	checkMulti(newcardname,inc_cats)


