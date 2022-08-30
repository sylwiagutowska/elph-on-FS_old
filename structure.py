Ha_to_ev=13.605662285137*2
PRECIS=6
import xml.etree.ElementTree as ET
import numpy as np
from operator import itemgetter
import glob

class structure():
 def __init__(self):
  self.allk=[]
  self.NONEQ=[]
  self.WK=[] #weights of kpoints
  self.pm=[0,1,-1]
  self.at=[]
  self.e=[]
  self.SYMM=[]
  self.no_of_kpoints=[]
  self.prefix=''
  self.prefixdyn=''
  self.at_pos=[]
  self.at_masses=[]
  self.at_names=[]
  self.irt=[] #we rotate atom by sym and check which atom he becomes
  self.rtau=[] #pos of atoms rotated by sym
  self.tmp_dir='../tmp_dir'
 
 def read_structure(self):
  print(' read structure...')
  self.prefix=glob.glob(self.tmp_dir+'/*.xml')[0].replace('/','.').split('.')[-2]
  self.prefixdyn=glob.glob('../*dyn0')[0].split('.dyn')[0]
  print(self.prefixdyn)
  tree = ET.parse(self.tmp_dir+'/'+self.prefix+'.xml')
  root = tree.getroot()

  for i in root.findall('output/band_structure/starting_k_points/monkhorst_pack'):
   self.no_of_kpoints=[ int(i.get('nk1')),int(i.get('nk2')),int(i.get('nk3'))]
  print('Monkhorst grid of '+str( self.no_of_kpoints[0])+' '+str( self.no_of_kpoints[1])+' '+str( self.no_of_kpoints[2])+' points')

  mmm=0
  for i in root.findall('output/band_structure/ks_energies'):
    for actor in i.findall('k_point'):
     self.NONEQ.append([round(float(m),PRECIS) for m in  actor.text.split()])
     self.NONEQ[-1].append(mmm)
     mmm=mmm+1
  print ('No of nonequivalent kpoints= '+str(len(self.NONEQ)))

  #reciprocal lattice vectors
  for i in root.findall('output/basis_set/reciprocal_lattice'):
   self.e.append([round(float(m),PRECIS) for m in i.find('b1').text.split()])
   self.e.append([round(float(m),PRECIS) for m in i.find('b2').text.split()])
   self.e.append([round(float(m),PRECIS) for m in i.find('b3').text.split()])
  self.e=np.array(self.e)
  print (self.e)

  # lattice vectors
  for i in root.findall('output/atomic_structure/cell'):
   self.at.append([round(float(m),PRECIS) for m in i.find('a1').text.split()])
   self.at.append([round(float(m),PRECIS) for m in i.find('a2').text.split()])
   self.at.append([round(float(m),PRECIS) for m in i.find('a3').text.split()])
  self.at=np.array(self.at)
  print (self.at)


  # atomic positions
  for i in root.findall('output/atomic_structure/atomic_positions/atom'):
   self.at_pos.append([round(float(m),PRECIS) for m in i.text.split()])
   self.at_names.append(i.get('name'))
  self.at_pos=np.array(self.at_pos)
  self.at_pos_crystal=[np.round(np.dot(np.linalg.inv(np.transpose(self.at)),i),PRECIS) for i in self.at_pos]
  for i in self.at_pos_crystal: 
    for m in range(3): 
     while i[m]>=1: i[m]=i[m]-1.
     while i[m]<0: i[m]=i[m]+1.
  print (self.at_pos) #,self.at_names)

  masses,names=[],[]
  #atomic masses
  for i in root.findall('input/atomic_species/species'):
   masses.append(round(float(i.find('mass').text),PRECIS)) #of spieces, not of positions
   names.append(i.get('name'))
  #print(masses,names)
  for i in self.at_names: # all positions
   for nj,j in enumerate(names): #all species
    if i==j: self.at_masses.append(masses[nj])

  #symmetry operations
  self.SYMM_crystal=[]
  self.SYMM=[]
  self.FRACTIONAL_TRANSLATION_crystal=[]
  self.FRACTIONAL_TRANSLATION_r=[]
  self.FRACTIONAL_TRANSLATION_k=[]
 # self.SYMM2=[]
  einv=np.linalg.inv(self.e)
 # atinv=np.linalg.inv(self.at)
  for neighbor in root.iter('rotation'):
     tmp=neighbor.text.split()
     tmp2=np.array([ [ float(m) for m in tmp[0:3]], [float(m) for m in tmp[3:6]], [float(m) for m in tmp[6:9]]])
     self.SYMM.append(np.transpose(np.dot(einv,(np.dot(tmp2,self.e)))))
     self.SYMM_crystal.append(tmp2)
  for neighbor in root.iter('fractional_translation'):
     tmp=neighbor.text.split()
     tmp2=np.round(np.array( [ float(m) for m in tmp[0:3]]),PRECIS)
     self.FRACTIONAL_TRANSLATION_crystal.append(tmp2)
     self.FRACTIONAL_TRANSLATION_r.append([round(sum([tmp2[m]*self.at[m][m2] for m in range(3)]),PRECIS) for m2 in range(3)]) #in cart --- depends wether in reciprocal or real space
     self.FRACTIONAL_TRANSLATION_k.append([round(sum([tmp2[m]* self.e[m][m2] for m in range(3)]),PRECIS) for m2 in range(3)]) #in cart --- depends wether in reciprocal or real space
  print ('No of symm. op.='+str(len(self.SYMM)))




  self.e=np.array(self.e)
  print (self.e)


 def calc_irt(self):
  einv=np.round(np.linalg.inv(np.transpose(self.at)),PRECIS)
  self.irt=[[-1 for i in self.SYMM] for j in self.at_pos]
  self.rtau=[[[] for i in self.SYMM] for j in self.at_pos]
  self.rtau_cart=[[[] for i in self.SYMM] for j in self.at_pos]
  for ni,i in enumerate(self.at_pos_crystal):
   for nj,j in enumerate(self.SYMM_crystal):
    pos22=[round(sum([j[m1][m2]*i[m2] for m2 in range(3)])+self.FRACTIONAL_TRANSLATION_crystal[nj][m1],PRECIS) for m1 in range(3)] #in cryst
    for m in range(3): 
     while pos22[m]>=1: pos22[m]=pos22[m]-1.
     while pos22[m]<0: pos22[m]=pos22[m]+1.
    pos2=np.round(pos22,PRECIS-1) #in cryst
    self.rtau[ni][nj]=pos2
    self.rtau_cart[ni][nj]=np.array([round(sum([pos2[m]*self.at[m][m2] for m in range(3)]),PRECIS) for m2 in range(3)])
    for nk,k in enumerate(self.at_pos_crystal):  
      pos2=np.round(pos2,PRECIS-2)
      k2=np.round(k,PRECIS-2)
      if pos2[0]==k2[0] and pos2[1]==k2[1] and pos2[2]==k2[2]:
       self.irt[ni][nj]=nk
       break
    if self.irt[ni][nj]==-1: raise ValueError('WRONG SYMMETRY: IRT NOT FOUND',j,pos2,pos22,self.at_pos_crystal)
 

 def sorting(self,allk2):
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
  
#  print ' Sorting - Done!'
  return allk

 def remove_repeated_items(self,allk2):
  found=1
  while found:
   allk=[]
   for i in range(len(allk2)-1):
    x=allk2[i]
    y=allk2[i+1]
    if not((x[0]==y[0]) and (x[1]==y[1]) and (x[2]==y[2])):
     allk.append(x)
   if not((allk[-1][0]==allk2[-1][0]) and (allk[-1][1]==allk2[-1][1]) and (allk[-1][2]==allk2[-1][2])): allk.append(allk2[-1])
   if len(allk)==len(allk2): found=0
   else: allk2=list(allk)
  return allk
 

 def calc_weight_of_k(self):
  WK=[0 for i in range(len(self.NONEQ))]
  for i in range(len(self.NONEQ)):
   for j in self.allk:
    if j[3]==i:
     WK[i]+=1
  return WK

 def make_kgrid(self,q=-1):
  print(' make whole kgrid')
  allk2=[]
  einv=np.linalg.inv(np.transpose(self.e))
  for nq in self.NONEQ:
  # print nq
#   nq2=[ round(sum([nq[m]*self.e[m][m2] for m in range(3)]),PRECIS)\
#                      for m2 in range(3)] #transform  to cartesian coordinates
   for ns,sym in enumerate(self.SYMM):
    x=[sum([sym[m1][m2]*nq[m2] for m2 in range(3)]) for m1 in range(3)]
    for h1 in self.pm:
     for k1 in self.pm:
      for l1 in self.pm:
       k_point2=[round(kk,PRECIS) for kk in 
                 (x-(h1*self.e[0]+k1*self.e[1]+l1*self.e[2]))]
       k_point3=[round(self.no_of_kpoints[m2]*\
                round(sum([k_point2[m]*einv[m2][m] for m in range(3)]),PRECIS)\
                     ) for m2 in range(3)] #transform from cartesian to crystal coordinates, then we have a cube of points
       if k_point3[0]>=0 and k_point3[0]<self.no_of_kpoints[0] and \
         k_point3[1]>=0 and k_point3[1]<self.no_of_kpoints[1] and \
         k_point3[2]>=0 and k_point3[2]<self.no_of_kpoints[2]:
         allk2.append(k_point2)
         allk2[-1].append(nq[3])
 # print (len(allk2))
  allk2=self.sorting(allk2)
  self.allk=allk2

  #rearrange
  allk_in_crystal_coordinates=[ [ round(self.no_of_kpoints[m2]*round(sum([v[m]*einv[m2][m] for m in range(3)]),PRECIS)) for m2 in range(3)]+[v[3],nv] for nv,v in enumerate(self.allk)]
  allk_in_crystal_coordinates=self.sorting(allk_in_crystal_coordinates)
  self.allk_in_crystal_coordinates=self.remove_repeated_items(allk_in_crystal_coordinates)

  self.allk=[ self.allk[i[4]] for i in self.allk_in_crystal_coordinates]

  self.allk=[i[:4] for i in self.allk]
  self.allk_in_crystal_coordinates=[i[:4] for i in self.allk_in_crystal_coordinates]
  for ni,i in enumerate(self.allk):
   self.allk_in_crystal_coordinates[ni][3]=i[3]


  #calc weights of k
  self.WK=self.calc_weight_of_k()


  if len(self.allk)!=self.no_of_kpoints[0]*self.no_of_kpoints[1]*self.no_of_kpoints[2]: 
   for i in range(self.no_of_kpoints[0]):  
    for j in range(self.no_of_kpoints[1]):
     for k in range(self.no_of_kpoints[2]):
      found=0
      for kk in self.allk_in_crystal_coordinates:
       if i==kk[0] and j==kk[1] and k==kk[2]: 
        found=1
        break
      if not found: print(q,': kpoint',i,j,k,'not found')
   raise ValueError(q,'wrong no of kpoints',len(self.allk),'!=',self.no_of_kpoints[0]*self.no_of_kpoints[1]*self.no_of_kpoints[2])

  '''
  h=open('kpoints.dat','w')
  for i in self.allk:
   for j in i:
    h.write(str(j)+' ')
   h.write('\n')
  h.close()

  h=open('kpoints0.dat','w')
  for i in allk_in_crystal_coordinates:
   for j in i:
    h.write(str(j)+' ')
   h.write('\n')
  h.close()

  h=open('kpoints_noneq.dat','w')
  for i in self.NONEQ:
   for j in i:
    h.write(str(j)+' ')
   h.write('\n')
  h.close()
  '''

 def check_symm(self,q,basic_kpoints):
  allk2=[]
 # SYMM2_crystal=[]
 # SYMM2=[]
  founded=[]
  einv=np.linalg.inv(np.transpose(self.e))
  k_pointp2=[round(self.no_of_kpoints[m2]*\
                round(sum([q[m]*einv[m2][m] for m in range(3)]),PRECIS)\
                     ) for m2 in range(3)] #transform from cartesian to crystal coordinates, then we have a cube of points  
 # print('cghecksymm',k_pointp2)
  for si,sym in enumerate(self.SYMM):
     found=0
     x=[sum([sym[m1][m2]*q[m2] for m2 in range(3)]) for m1 in range(3)]   #in cart
     for h1 in self.pm:
      for k1 in self.pm:
       for l1 in self.pm:
        k_point2=[round(kk,PRECIS) for kk in 
                 (x-(h1*self.e[0]+k1*self.e[1]+l1*self.e[2]))]
        k_point3=[round(self.no_of_kpoints[m2]*\
                round(sum([k_point2[m]*einv[m2][m] for m in range(3)]),PRECIS)\
                     ) for m2 in range(3)] #transform from cartesian to crystal coordinates, then we have a cube of points
        if k_point3[0]==k_pointp2[0] and  k_point3[1]==k_pointp2[1] and  k_point3[2]==k_pointp2[2]:
           found=1
           break
       if found: break
      if found: break
     if found: 
      founded.append(si)
#      SYMM2.append(sym)
#      SYMM2_crystal.append(self.SYMM_crystal[si])
  print(q,':',len(founded),'operations of symmetry')
  self.SYMM=[self.SYMM[si] for si in founded]
  self.SYMM_crystal=[self.SYMM_crystal[si] for si in founded]
  self.FRACTIONAL_TRANSLATION_crystal=[self.FRACTIONAL_TRANSLATION_crystal[si] for si in founded]
  self.FRACTIONAL_TRANSLATION_k=[self.FRACTIONAL_TRANSLATION_k[si] for si in founded]
  self.FRACTIONAL_TRANSLATION_r=[self.FRACTIONAL_TRANSLATION_r[si] for si in founded]

 def find_k_plus_q(self,k,allk,allk_cryst,q):
  kpq_no=-2
  einv=np.linalg.inv(np.transpose(self.e))
  k2=[round(sum([k[m]*einv[m2][m] for m in range(3)]),PRECIS) for m2 in range(3)]
  q2=[round(sum([q[m]*einv[m2][m] for m in range(3)]),PRECIS) for m2 in range(3)]
  kpq0=[k2[i]+q2[i] for i in range(3)] #in cryst coords
  for i in range(3): 
   while kpq[i]>=1: kpq[i]=kpq[i]-1
   while kpq[i]<0: kpq[i]=kpq[i]+1

  kpq=[round(sum([kpq0[m]*self.e[m][m2] for m in range(3)]),PRECIS) for m2 in range(3)] #in cart
  found=0
  for sym in self.SYMM:
     if found==1: break
     x=[sum([sym[m1][m2]*kpq[m2] for m2 in range(3)]) for m1 in range(3)] #in cart coord
     kpq2=[round(self.no_of_kpoints[m2]*\
                 round(sum([x[m]*einv[m2][m] for m in range(3)]),PRECIS)\
                     ) for m2 in range(3)] #transform from cartesian to crystal coordinates, then we have a cube of points
  for i in range(3): 
   while kpq2[i]>=self.no_of_kpoints[i]: kpq2[i]=kpq2[i]-self.no_of_kpoints[i]
   while kpq2[i]<0: kpq2[i]=kpq2[i]+self.no_of_kpoints[i]
  for ki in allk_cryst:
      if found==1: break
      if kpq2[0]==ki[0] and kpq2[1]==ki[1] and kpq2[2]==ki[2]:
            kpq_no=ki[3]
            kpq_no_in_basic=ki[4]
            found=1
            break
  if kpq_no==-2: raise ValueError(q,k,kpq,' kpq not found')
  return [kpq,kpq_no,kpq_no_in_basic]

 def find_newkpoints_in_old_list(self,old_allk):
  for i in self.NONEQ: i.append(None)
  for nk,k in enumerate(self.allk):
   k.append(old_allk[nk][3]) # k=kx,ky,kz,no_of_noneq_in_new_grid,no_of_noneq_in_basic_grid
   self.allk_in_crystal_coordinates[nk].append(old_allk[nk][3])
   self.NONEQ[k[3]][4]=old_allk[nk][3]

