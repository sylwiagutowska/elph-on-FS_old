import xml.etree.ElementTree as ET
Ha_to_ev=13.605662285137*2

class el_structure():
 def __init__(self,structure):
  self.ENE_fs=[]
  self.ENE=[]
  self.ef=0.
  self.bands_num=[]
  self.minband=0
  self.maxband=0
  self.prefix=structure.prefix

 def read_el_structure(self):
  tree = ET.parse('tmp_dir/'+self.prefix+'.xml')
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
  self.minband,self.maxband=self.bands_num[0],self.bands_num[-1]

  self.ENE=ENE2
  self.ENE_fs=ENE2[self.minband:self.maxband+1] #np.transpose(np.array(ENE2)) #ENE=[] #ENE[i][j] , i - no of band, j-no of kpoint


