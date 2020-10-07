dir='tmp_dir/_ph0/ir.phsave/'
import xml.etree.ElementTree as ET
import numpy as np
from operator import itemgetter

PRECIS=6
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
ENE=[] #ENE[i][j] , i -no of kpoint, j-no of band
NONEQ=[] #NONEQ[0:3]=kpoint, NONEQ[3]=no of kpoint

for i in root.findall('output/band_structure/fermi_energy'):
 ef=float(i.text.split()[0])
print(ef)

for i in root.findall('output/band_structure/starting_k_points/monkhorst_pack'):
 no_of_kpoints=[ int(i.get('nk1')),int(i.get('nk2')),int(i.get('nk3'))]
print( no_of_kpoints)


mmm=0
for i in root.findall('output/band_structure/ks_energies'):
    for actor in i.findall('eigenvalues'):
     ENE.append([float(m)-ef for m in  actor.text.split()])
    for actor in i.findall('k_point'):
     NONEQ.append([round(float(m),PRECIS) for m in  actor.text.split()])
     NONEQ[-1].append(mmm)
     mmm=mmm+1

maks=max([ len(i) for i in ENE])
ENE2= [ [] for i in range(maks)]
below_ENE=[ 0 for i in range(maks)]
top_ENE=[ 0 for i in range(maks)]

#choose only complete bands
for i in ENE:
 for (numj,j) in enumerate(i):
  ENE2[numj].append(j)
  if j<0: below_ENE[numj]=1
  elif j>0: top_ENE[numj]=1
bands_num=[]
#print('\nList of bands, which cross the Fermi level and its energy ranges:')
for i in range(maks):
 if below_ENE[i] and top_ENE[i]: 
  bands_num.append(i)
minband,maxband=bands_num[0],bands_num[-1]

ENE=ENE2[minband:maxband+1] #np.transpose(np.array(ENE2)) #ENE=[] #ENE[i][j] , i - no of band, j-no of kpoint

#reciprocal lattice vectors
e=[]
for i in root.findall('output/basis_set/reciprocal_lattice'):
 e.append([round(float(m),PRECIS) for m in i.find('b1').text.split()])
 e.append([round(float(m),PRECIS) for m in i.find('b2').text.split()])
 e.append([round(float(m),PRECIS) for m in i.find('b3').text.split()])
e=np.array(e)
print e

#symmetry operations
SYMM=[np.array([[1,0,0],[0,1,0],[0,0,1]])]
for neighbor in root.iter('rotation'):
     tmp=neighbor.text.split()
     tmp2=np.array([ [ float(m) for m in tmp[0:3]], [float(m) for m in tmp[3:6]], [float(m) for m in tmp[6:9]]])
     SYMM.append(np.transpose(np.dot((np.linalg.inv(e)),(np.dot(tmp2,e)))))
print len(SYMM)
### 



###make whole kgrid

allk=[]
mmm=0
pm=[-1,0.,1]
A_vectors=[1.1*np.cross(e[0],e[1]),np.cross(e[2],e[0]),np.cross(e[1],e[2]),\
           -np.cross(e[0],e[1]),-np.cross(e[2],e[0]),-np.cross(e[1],e[2])]
B_vectors=[e[2]+0.5*(e[0]+e[1]),e[1]+0.5*(e[0]+e[2]),e[0]+0.5*(e[1]+e[2]),\
0.5*(e[0]+e[1]),0.5*(e[0]+e[2]),0.5*(e[1]+e[2]) ]
A_vectors=[ [round(kk,PRECIS) for kk in m] for m in A_vectors]
B_vectors=[ [round(kk,PRECIS) for kk in m] for m in B_vectors]

#NONEQ=[ [ round(sum([v[m]*np.transpose((e))[m2][m] for m in range(3)]),4) for m2 in range(3)]+[v[3]] for v in NONEQ]
NONEQ2=[]
einv=np.linalg.inv(np.transpose(e))
for nq in NONEQ:
 # print nq
  for sym in SYMM:
   x=[sum([sym[m1][m2]*nq[m2] for m2 in range(3)]) for m1 in range(3)]
   for h1 in pm:
    for k1 in pm:
     for l1 in pm:
      k_point2=[round(kk,PRECIS) for kk in (x-(h1*e[0]+k1*e[1]+l1*e[2]))]
      k_point3=[round(no_of_kpoints[m2]*\
               round(sum([k_point2[m]*einv[m2][m] for m in range(3)]),PRECIS)\
                     ) for m2 in range(3)] #transform from cartesian to crystal coordinates, then we have a cube of points
      if k_point3[0]>=0 and k_point3[0]<=no_of_kpoints[0] and \
         k_point3[1]>=0 and k_point3[1]<=no_of_kpoints[1] and \
         k_point3[2]>=0 and k_point3[2]<=no_of_kpoints[2]:
         allk.append(k_point2)
         allk[-1].append(nq[3])
      '''
      k_point2=x-(h1*e[0]+k1*e[1]+l1*e[2])
      sign0=0
      for numa,a in enumerate(A_vectors):
        ra=k_point2-B_vectors[numa]
        scalar_product=round(sum([a[i]*ra[i] for i in range(3)]),PRECIS)
        if scalar_product>0: #  or (numa<3 and scalar_product==0): 
         sign0=1
         break
      if sign0==0:
        allk.append([round(kk,PRECIS) for kk in k_point2])
        allk[-1].append(nq[3])
      '''
print len(allk)
allk2=sorting(allk)

allk=[] 
for i in range(len(allk2)-1):
  x=allk2[i]
  y=allk2[i+1]
  if not((x[0]==y[0]) and (x[1]==y[1]) and (x[2]==y[2])):
   allk.append(x)
print len(allk)

h=open('kpoints.dat','w')
for i in allk:
 for j in i:
  h.write(str(j)+' ')
 h.write('\n')
h.close()

#rearrange
allk2=[ [ round(no_of_kpoints[0]*round(sum([v[m]*einv[m2][m] for m in range(3)]),PRECIS)) for m2 in range(3)]+[nv] for nv,v in enumerate(allk)]
allk2=sorting(allk2)
allk=[ allk[i[3]] for i in allk2 if i[0]<no_of_kpoints[0]  and i[1]<no_of_kpoints[0] and i[2]<no_of_kpoints[0]]

#allk=sorting(allk)
print(len(allk))

h=open('kpoints0.dat','w')
for i in allk2:
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
for file in range(2,3):
 for mode in range(1,len(NONDEG[file-1])+1):
  KPOINTS,ELPH,COLORS=[],[],[]
  tree = ET.parse(dir+'elph.'+str(file)+'.'+str(mode)+'.xml')
  root = tree.getroot()
  for country in root.iter('PARTIAL_EL_PHON'):
   nkp=int(country.find('NUMBER_OF_K').text)
   nbnd_el=int(country.find('NUMBER_OF_BANDS').text)
   for k in range(1,nkp+1):
    for town in country.iter('K_POINT.'+str(k)):
     KPOINTS.append([ round(float(m),PRECIS) for m in town.find('COORDINATES_XK').text.split() ])
  #   print town.find('PARTIAL_ELPH').text.split('\n')
     elph_k=([  (float(m.replace(',',' ').split()[0])**2+float(m.replace(',',' ').split()[1])**2) for m in town.find('PARTIAL_ELPH').text.split('\n') if len(m.split())>0  ])
     ELPH.append([ sum(elph_k[nbnd_el*i:nbnd_el*(i+1)]) for i in bands_num]) #choose only bands which cross EF and sum over j in pairs <i,j>

  
#  print KPOINTS
  for num_noneqk,noneqk in enumerate(NONEQ):
   colored=0
   for numk,k in enumerate(KPOINTS):
    if (k[0]==noneqk[0]) and  (k[1]==noneqk[1]) and  (k[2]==noneqk[2]):
     COLORS.append(ELPH[numk])
     colored=1
     break
   if colored==0:
    COLORS.append([0 for i in bands_num])
  COLORS=np.transpose(np.array(COLORS)) #COLORS[nbnd][nkp]
  print nbnd_el,bands_num
  print len(COLORS),[len(i) for i in COLORS]
  print len(ENE),[len(i) for i in ENE]
  print len(allk)
  h=open('elph'+str(mode)+'.frmsf','w')
  h.write(str(no_of_kpoints[0])+' '+str(no_of_kpoints[1])+' ' +str(no_of_kpoints[2])+'\n')
  h.write('1\n'+str(len(bands_num))+'\n')
  for i in e:
   for j in i:
    h.write(str(j)+' ')
   h.write('\n')
  for bnd in ENE:
   for k in allk:
    h.write(str(bnd[k[3]])+'\n')
  for bnd in COLORS:
   for k in allk:
    h.write(str(bnd[k[3]])+'\n')
  h.close()



###
