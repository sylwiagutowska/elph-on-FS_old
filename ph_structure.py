import xml.etree.ElementTree as ET
PRECIS=6
class ph_structure():
 def __init__(self):
  self.FREQ=[] #[i][j] i: q-point; j-nbnd
  self.NONDEG=[]
  self.DEG=[] #DEG[q][nbnd] = no of band, with whom the nbnd is degenerated
  self.nat=0
  self.no_of_modes=0
  self.no_of_atoms=0
  self.ELPH_sum=[]
  self.elph_dir='tmp_dir/_ph0/ir.phsave/'
 def read_ph_structure(self):
  ###read frequencies and make DEG[q][nbnd] = no of band, with whom the nbnd is degenerated
  h=open('ir.dyn1')
  for i in range(3): tmp=h.readline()
  h.close()
  self.nat=int(tmp.split()[1])

  for file in range(1,30):
   h=open('ir.dyn'+str(file))
   tmp=h.readlines()
   h.close()
   self.FREQ.append([])
   for i in tmp:
    if 'freq' in i:
     self.FREQ[-1].append(float(i.split()[4]))
  self.no_of_modes=len(self.FREQ[-1])
  self.no_of_atoms=int(self.no_of_modes/3)
    
  for f in self.FREQ:
   self.NONDEG.append([])
   self.DEG.append([0 for i in range(3*self.nat)])
   for j in range(3*self.nat):
    is_deg=0
    for numf2,f2 in enumerate(self.NONDEG[-1]):
     if abs(f[j]-f2)<1e-3:
      self.DEG[-1][j]=numf2
      is_deg=1
      break
    if is_deg==0:
     self.NONDEG[-1].append(f[j])
     self.DEG[-1][j]=len(self.NONDEG[-1])-1

       
