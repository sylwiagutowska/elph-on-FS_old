dir='tmp_dir/_ph0/ir.phsave/'
import xml.etree.ElementTree as ET
import numpy as np
from operator import itemgetter

PRECIS=5
def sorting(allk2):
 Xall=[]
 allk2=sorted(allk2, key=itemgetter(0))
 i=0
 while i<len(allk2): 
  X=[]
  x=allk2[i]
  while i<len(allk2) and x[0]==allk2[i][0]:
   X.append(allk2[i])
   i=i+1
  if len(X)>1: X=sorted(X, key=itemgetter(1))
  Xall.append(X)

 Yall=[]
 for X in Xall:
  x=X[0]
  i=0
  while i<len(X): 
   Y=[]
   x=X[i]
   while i<len(X) and x[1]==X[i][1]:
    Y.append(X[i])
    i=i+1
   Y=sorted(Y, key=itemgetter(2))
   Yall.append(Y)

 allk=[]
 for i in Yall:
  for j in i:
   allk.append(j)

 print ' Sorting - Done!'

 return allk

###read energies and other electronic and structure info
Ha_to_ev=13.605662285137*2
#file=raw_input('Path and name of *.xml file, eq. Pb or tmp_dir/Pb: ')
tree = ET.parse('FS/tmp_dir/ir.xml')
root = tree.getroot()
ENE=[]
NONEQ=[]

for i in root.findall('output/band_structure/fermi_energy'):
 ef=float(i.text.split()[0])
print(ef)

for i in root.findall('output/band_structure/starting_k_points/monkhorst_pack'):
 no_of_kpoints=[ int(i.get('nk1')),int(i.get('nk2')),int(i.get('nk3'))]
print( no_of_kpoints)

for i in root.findall('output/band_structure/ks_energies'):
    for actor in i.findall('eigenvalues'):
     ENE.append([float(m) for m in  actor.text.split()])
    for actor in i.findall('k_point'):
     NONEQ.append([round(float(m),PRECIS) for m in  actor.text.split()])

maks=max([ len(i) for i in ENE])
ENE2= [ [] for i in range(maks)]
below_ENE=[ 0 for i in range(maks)]
top_ENE=[ 0 for i in range(maks)]

for i in ENE:
 for (numj,j) in enumerate(i):
  ENE2[numj].append(j)
  if j<ef: below_ENE[numj]=1
  elif j>ef: top_ENE[numj]=1
bands_num=[]
#print('\nList of bands, which cross the Fermi level and its energy ranges:')
for i in range(maks):
 if below_ENE[i] and top_ENE[i]: 
  bands_num.append(i)
minband,maxband=bands_num[0],bands_num[-1]

#reciprocal lattice vectors
e=[]
for i in root.findall('output/basis_set/reciprocal_lattice'):
 e.append([round(float(m),PRECIS) for m in i.find('b1').text.split()])
 e.append([round(float(m),PRECIS) for m in i.find('b2').text.split()])
 e.append([round(float(m),PRECIS) for m in i.find('b3').text.split()])
e=np.transpose(np.array(e))
print e

#symmetry operations
SYMM=[np.array([[1,0,0],[0,1,0],[0,0,1]])]
for neighbor in root.iter('rotation'):
     tmp=neighbor.text.split()
     tmp2=np.array([ [ float(m) for m in tmp[0:3]], [float(m) for m in tmp[3:6]], [float(m) for m in tmp[6:9]]])
     SYMM.append(np.dot((np.linalg.inv(e)),(np.dot((tmp2),e))))
print len(SYMM)
### 



###make whole kgrid

allk=[]
mmm=0
pm=[-1,0.,1,-2,2]



NONEQ=[ [ round(sum([v[m]*((e))[m2][m] for m in range(3)]),4) for m2 in range(3)] for v in NONEQ]
#NONEQ=[ [ (sum([v[m]*np.linalg.inv((e))[m2][m] for m in range(3)])) for m2 in range(3)] for v in NONEQ]
#NONEQ=[ [ round(sum([v[m]*np.transpose((e))[m2][m] for m in range(3)]),4) for m2 in range(3)] for v in NONEQ]
NONEQ=[ [ int(no_of_kpoints[0]/2*v[m]) for m in range(3)] for v in NONEQ]
e2=no_of_kpoints[0]*np.array([[1,0,0],[0,1,0],[0,0,1]])
NONEQ=sorting(NONEQ)
'''
NONEQ2=[]
#trans
print len(NONEQ)

for nq in NONEQ:
   for h1 in pm:
    for k1 in pm:
     for l1 in pm:
      k_point2=[int(m) for m in nq+(h1*e2[0]+k1*e2[1]+l1*e2[2]) ]
      sign0=0
      if k_point2[0]>=0 and k_point2[0]<=no_of_kpoints[0] and k_point2[1]>=0 and k_point2[1]<=no_of_kpoints[0] and k_point2[2]>=0 and k_point2[2]<=no_of_kpoints[0]:
       NONEQ2.append(k_point2)

NONEQ2=sorting(NONEQ2)

NONEQ=[] 
for i in range(len(NONEQ2)-1):
  x=NONEQ2[i]
  y=NONEQ2[i+1]
  if not(x[0]==y[0] and x[1]==y[1] and x[2]==y[2]):
   NONEQ.append(x)
print len(NONEQ)
'''

for nq in NONEQ:
 # print nq
  for sym in SYMM:
   x=np.array([int(sum([sym[m2][m1]*nq[m2] for m2 in range(3)])) for m1 in range(3)])
   #if x[0]<maxe[0] and x[1]<maxe[1] and x[2]<maxe[2] and x[0]>-maxe[0] and x[1]>-maxe[1] and x[2]>-maxe[2]:
   for h1 in pm:
    for k1 in pm:
     for l1 in pm:
      k_point2=[int(m) for m in x-(h1*e2[0]+k1*e2[1]+l1*e2[2]) ]
      if k_point2[0]>=-no_of_kpoints[0]/2 and k_point2[0]<no_of_kpoints[0]/2 and k_point2[1]>=-no_of_kpoints[0]/2 and k_point2[1]<no_of_kpoints[0]/2 and k_point2[2]>=-no_of_kpoints[0]/2 and k_point2[2]<no_of_kpoints[0]/2:
       allk.append(k_point2)

print len(allk)
allk2=sorting(allk)

allk=[] 
for i in range(len(allk2)-1):
  x=allk2[i]
  y=allk2[i+1]
  if not(x[0]==y[0] and x[1]==y[1] and x[2]==y[2]):
   allk.append(x)
print len(allk)



h=open('kpoints.dat','w')
for i in allk:
 for j in i:
  h.write(str(j)+' ')
 h.write('\n')
h.close()

allk=[ [ round(sum([v[m]*np.linalg.inv(np.transpose(e))[m2][m] for m in range(3)]),4) for m2 in range(3)] for v in allk]
allk=sorting(allk)


h=open('kpoints0.dat','w')
for i in allk:
 for j in i:
  h.write(str(j)+' ')
 h.write('\n')
h.close()

h=open('kpoints_noneq.dat','w')
for i in NONEQ:
 for j in i:
  h.write(str(j)+' ')
 h.write('\n')
h.close()

###

###read frequencies and make DEG[q][nbnd] = no of band, with whom the nbnd is degenerated
h=open('ir.dyn1')
for i in range(3): tmp=h.readline()
h.close()
nat=int(tmp.split()[1])

FREQ=[] #[i][j] i: q-point; j-nbnd
for file in range(1,30):
 h=open('ir.dyn'+str(file))
 tmp=h.readlines()
 h.close()
 FREQ.append([])
 for i in tmp:
  if 'freq' in i:
   FREQ[-1].append(float(i.split()[4]))

NONDEG=[]
DEG=[]
for f in FREQ:
 NONDEG.append([])
 DEG.append([0 for i in range(3*nat)])
 for j in range(3*nat):
  is_deg=0
  for numf2,f2 in enumerate(NONDEG[-1]):
   if abs(f[j]-f2)<1e-3:
    DEG[-1][j]=numf2
    is_deg=1
    break
  if is_deg==0:
   NONDEG[-1].append(f[j])
   DEG[-1][j]=len(NONDEG[-1])-1

####

###read elph matrix
for file in range(1,2):
 for j in range(1,len(NONDEG[file-1])+1):
  KPOINTS,ELPH,COLORS=[],[],[]
  tree = ET.parse(dir+'elph.'+str(file)+'.'+str(j)+'.xml')
  root = tree.getroot()
  for country in root.iter('PARTIAL_EL_PHON'):
   nkp=int(country.find('NUMBER_OF_K').text)
   nbnd_el=int(country.find('NUMBER_OF_BANDS').text)
   nbnd_el=int(country.find('NUMBER_OF_BANDS').text)
   for k in range(1,nkp+1):
    for town in country.iter('K_POINT.'+str(k)):
     KPOINTS.append([ float(m) for m in town.find('COORDINATES_XK').text.split() ])
  #   print town.find('PARTIAL_ELPH').text.split('\n')
     elph_k=([ [ float(m.replace(',',' ').split()[0])**2+float(m.replace(',',' ').split()[1])**2 for m in town.find('PARTIAL_ELPH').text.split('\n') if len(m.split())>0 ] ])
     ELPH.append([ sum(elph_k[nbnd_el*i:nbnd_el*(i+1)]) for i in bands_num]) #choose only bands which cross EF and sum over j in pairs <i,j>


  for num_noneqk,noneqk in enumerate(NONEQ):
   colored=0
   for numk,k in enumerate(KPOINTS):
    if (k[0]-noneqk[0])<1e-3 and  (k[1]-noneqk[1])<1e-3 and  (k[2]-noneqk[2])<1e-3:
     COLORS.append(ELPH[numk])
     colored=1
     break
   if colored==0:
    COLORS.append([0 for i in ELPH[numk]])


 print len(COLORS)

###
