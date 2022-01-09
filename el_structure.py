import xml.etree.ElementTree as ET
from numpy import array,reshape
Ha_to_ev=13.605662285137*2

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
  self.minband,self.maxband=self.bands_num[0],self.bands_num[-1]
  self.fermi_nbnd_el=len(self.bands_num)

  self.ENE=ENE2
  self.ENE_fs=ENE2[self.minband:self.maxband+1] #np.transpose(np.array(ENE2)) #ENE=[] #ENE[i][j] , i - no of band, j-no of kpoint


 def read_extra_el_structure(self):
  ##read dense grid energt
  print(' read electronic structure on dense grid')
  data=[]
  ENE_extra=[]
  h=open(self.tmp_dir+'/'+self.prefix+'.a2Fsave')
  tmp=[i.split() for i in h.readlines()]
  h.close()
  [nbnd_extra,nkp_extra]=[int(m) for m in tmp[0]]
  for i in tmp[1:]:
   if len(i)==3: 
    self.nk_extra=[int(m) for m in i]
    break
   for m in i:
    data.append(float(m))
  print(len(data))
  self.ENE_extra=array(data[:nkp_extra*nbnd_extra]).reshape(nbnd_extra,nkp_extra)
  self.NONEQ_extra=array(data[nkp_extra*nbnd_extra:nkp_extra*nbnd_extra+nkp_extra*3]).reshape(nkp_extra,3)
  self.WK_extra=array(data[nkp_extra*nbnd_extra+nkp_extra*3:])
  self.NONEQ_extra=[list(m)+[nm] for nm,m in enumerate(self.NONEQ_extra)]

   #ENE_extra.append(np.array([float(m) for m in i]).reshape(nbnd_
#from elphon:
#     READ(iuna2Fsave,*) ibnd, nksfit
#     READ(iuna2Fsave,*) etfit
#     READ(iuna2Fsave,*) ((xkfit(i,ik), i=1,3), ik=1,nksfit)
#     READ(iuna2Fsave,*) wkfit
#     READ(iuna2Fsave,*) nk1fit, nk2fit, nk3fit
