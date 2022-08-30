import xml.etree.ElementTree as ET
from numpy import exp
Ha_to_ev=13.605662285137*2
def w0gauss(x):
  degauss=0.1
  x2=x/degauss
  sqrtpm1= 1. / 1.77245385090551602729
  # cold smearing  (Marzari-Vanderbilt-DeVita-Payne)
  arg = min (200., (x2 - 1.0 /  (2.0)**0.5 ) **2.)
  return sqrtpm1 * exp ( - arg) * (2.0 - ( 2.0)**0.5 * x2)/degauss

class el_structure():
 def __init__(self,structure):
  self.ENE_fs=[]
  self.ENE=[]
  self.ef=0.
  self.bands_num=[]
  self.fermi_nbnd_el=0
  self.minband=0
  self.maxband=0
  self.prefix=structure.prefix
  self.tmp_dir=structure.tmp_dir

 def read_el_structure(self):
  print(' read electronic structure')
  tree = ET.parse(self.tmp_dir+'/'+self.prefix+'.xml')
  root = tree.getroot()

  for i in root.findall('output/band_structure/fermi_energy'):
   self.ef=float(i.text.split()[0])
  print(self.ef)

  for i in root.findall('output/band_structure/ks_energies'):
    for actor in i.findall('eigenvalues'):
     self.ENE.append([float(m)-self.ef for m in  actor.text.split()])
  
  maks=max([ len(i) for i in self.ENE])
  ENE2= [ [] for i in range(maks)]
  below_ENE=[ 0 for i in range(maks)]
  top_ENE=[ 0 for i in range(maks)]

#choose only complete bands
  for i in self.ENE:
   for (numj,j) in enumerate(i):
    ENE2[numj].append(j)
    if j<0: below_ENE[numj]=1
    elif j>0: top_ENE[numj]=1
#print('\nList of bands, which cross the Fermi level and its energy ranges:')
  for i in range(maks):
   if below_ENE[i] and top_ENE[i]: 
    self.bands_num.append(i)
  self.bands_num=self.bands_num[::2]
  self.minband,self.maxband=self.bands_num[0],self.bands_num[-1]
  self.fermi_nbnd_el=len(self.bands_num)

  self.ENE=ENE2
  self.ENE_fs=ENE2[self.minband:self.maxband+1:2] #np.transpose(np.array(ENE2)) #ENE=[] #ENE[i][j] , i - no of band, j-no of kpoint


