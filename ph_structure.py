import xml.etree.ElementTree as ET
PRECIS=6
class ph_structure():
 def __init__(self,structure):
  self.FREQ=[] #[i][j] i: q-point; j-nbnd
  self.NONDEG=[]
  self.DEG=[] #DEG[q][nbnd] = no of band, with whom the nbnd is degenerated
  self.nat=0
  self.Q=[]
  self.no_of_modes=0
  self.DYN=[]
  self.multiplicity_of_qs=[]
#  self.SYMMQ=[]
  self.PATT=[]
  self.prefix=structure.prefix
  self.elph_dir=structure.tmp_dir+'/_ph0/'+self.prefix+'.phsave/'  

 def read_dyn_of_q(self,tmp):
  for ni,i in enumerate(tmp):
    if 'q = ' in i:
     tmp_dyn=[[float(x) for x in line.split()] for line in tmp[ni+2:ni+2+self.nat*self.nat*4]]
     tmp_dyn=[[[complex(ln[0],ln[1]), complex(ln[2],ln[3]), complex(ln[4],ln[5])] \
             for ln in tmp_dyn[m+1:m+4]] \
             for m in range(len(tmp_dyn)) if m%4==0]
     dyn=[ [0 for x in range(self.no_of_modes)] for y in range(self.no_of_modes)]
     for ii in range(3):
      for na in range(self.nat):
       mu=3*na+ii
       for jj in range(3):
        for nb in range(self.nat):
         nu=3*nb+jj
         dyn[mu][nu]=tmp_dyn[na*self.nat*3+nb*3][ii][jj] #/sqrt(mass_a*mass_b)/amu_ry
     break #if 'Diagonalizing' in i: break  
  self.DYN.append(dyn)

 def read_freq_of_q(self,tmp):
   self.FREQ.append([])
   for i in tmp:
    if 'freq' in i:
     self.FREQ[-1].append(float(i.split()[4]))
   self.no_of_modes=len(self.FREQ[-1])

 def read_q(self,tmp):
   for ni,i in enumerate(tmp):
    if 'q = ' in i:
     self.Q.append([round(float(m),PRECIS) for m in i.split()[3:6] ])
     break 

 def check_degeneration_of_modes(self):
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

 def read_multiplicity_of_q(self,dynfile):
  multiplicity=0
  for i in dynfile:
   if 'Dynamical  Matrix in cartesian axes' in i:
    multiplicity=multiplicity+1
  self.multiplicity_of_qs.append( multiplicity)

 def read_ph_structure(self):
  ###read frequencies and make DEG[q][nbnd] = no of band, with whom the nbnd is degenerated
  print(' read phonon structure')
  h=open(self.prefix+'.dyn1')
  for i in range(3): tmp=h.readline()
  h.close()
  self.nat=int(tmp.split()[1])
  for file in range(1,100):
   try: h=open(self.prefix+'.dyn'+str(file))
   except: break
   tmp=h.readlines()
   h.close()
   self.read_freq_of_q(tmp)
   self.read_dyn_of_q(tmp)
   self.read_q(tmp)
   self.read_multiplicity_of_q(tmp)
  self.check_degeneration_of_modes()



 '''
 def check_symm_of_q(self,structure):
  self.SYMMQ=[ ]
  self.pm=[0,1,-1]
  for nq in self.Q:
   self.SYMMQ.append([])
   for sym in structure.SYMM:
    found=0
    x=[sum([sym[m1][m2]*nq[m2] for m2 in range(3)]) for m1 in range(3)]
    for h1 in self.pm:
     for k1 in self.pm:
      for l1 in self.pm:
       q2=[round(kk,PRECIS) for kk in 
                 (x-(h1*structure.e[0]+k1*structure.e[1]+l1*structure.e[2]))]
       if (nq[0]==q2[0] and nq[1]==q2[1] and nq[2]==q2[2]):
         found=1
         break
      if found==1: break
     if found==1: 
      self.SYMMQ[-1].append(sym)
      break
#  for i in self.SYMMQ: print len(i)
 '''       

 def read_patterns(self):
  for q in range(len(self.Q)):
   self.PATT.append([])
   tree = ET.parse(self.elph_dir+'/patterns.'+str(q+1)+'.xml')
   root = tree.getroot()
   for i in range(len(self.NONDEG[q])):
     rep=root.find('IRREPS_INFO/REPRESENTION.'+str(i+1))
     npert=int(rep.find('NUMBER_OF_PERTURBATIONS').text)
     for j in range(npert):
      pat=[ (lambda m: complex(float(m[0]),float(m[1])))(n.replace(',',' ').split()) for n in  rep.find('PERTURBATION.'+str(j+1)+'/DISPLACEMENT_PATTERN').text.split('\n')[1:-1]]
      self.PATT[-1].append(pat)

     
