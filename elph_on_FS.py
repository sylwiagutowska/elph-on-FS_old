dir='tmp_dir/_ph0/ir.phsave/'
import xml.etree.ElementTree as ET

###read energies
Ha_to_ev=13.605662285137*2
#file=raw_input('Path and name of *.xml file, eq. Pb or tmp_dir/Pb: ')
tree = ET.parse('FS/tmp_dir/ir.xml')
root = tree.getroot()
ENE=[]
NONEQ=[]
#for child in root:
#     print child.tag, child.attrib
for i in root.findall('output/band_structure/fermi_energy'):
 ef=float(i.text.split()[0])
print(ef)

for i in root.findall('output/band_structure/ks_energies'):
    for actor in i.findall('eigenvalues'):
     ENE.append([float(m) for m in  actor.text.split()])
    for actor in i.findall('k_point'):
     NONEQ.append([float(m) for m in  actor.text.split()])

print (len(ENE))
print (len(NONEQ))

maks=max([ len(i) for i in ENE])
ENE2= [ [] for i in range(maks)]
below_ENE=[ 0 for i in range(maks)]
top_ENE=[ 0 for i in range(maks)]

for i in ENE:
 for (numj,j) in enumerate(i):
  ENE2[numj].append(j)
  if j<ef: below_ENE[numj]=1
  elif j>ef: top_ENE[numj]=1

#print('\nFermi level='+str(ef*Ha_to_ev))

#print('\nList of bands and its energy ranges:')
#for i in range(maks):
# print i+1, min(ENE2[i])*Ha_to_ev, max(ENE2[i])*Ha_to_ev

bands_num=[]
#print('\nList of bands, which cross the Fermi level and its energy ranges:')
for i in range(maks):
 if below_ENE[i] and top_ENE[i]: 
  bands_num.append(i)
minband,maxband=bands_num[0],bands_num[-1]
#  print i+1, min(ENE2[i])*Ha_to_ev, max(ENE2[i])*Ha_to_ev
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
