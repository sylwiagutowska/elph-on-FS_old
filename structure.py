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
  self.tmp_dir='../tmp_dir'

 def read_structure(self):
  print(' read structure...')
  self.prefix=glob.glob(self.tmp_dir+'/*.xml')[0].replace('/','.').split('.')[-2]
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

  #symmetry operations
  self.SYMM=[np.array([[1,0,0],[0,1,0],[0,0,1]])]
  einv=np.linalg.inv(self.e)
  for neighbor in root.iter('rotation'):
     tmp=neighbor.text.split()
     tmp2=np.array([ [ float(m) for m in tmp[0:3]], [float(m) for m in tmp[3:6]], [float(m) for m in tmp[6:9]]])
     self.SYMM.append(np.transpose(np.dot(einv,(np.dot(tmp2,self.e)))))
  print ('No of symm. op.='+str(len(self.SYMM)))

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
  allk=[]
  for i in range(len(allk2)-1):
    x=allk2[i]
    y=allk2[i+1]
    if not((x[0]==y[0]) and (x[1]==y[1]) and (x[2]==y[2])):
     allk.append(x)
  print (len(allk))
  return allk
 
 ''' doesnot work 
 def check_symm(self):
  SYMM2=[]
  einv=np.linalg.inv(np.transpose(self.e))
  for sym in self.SYMM:
   found=0
   for nq in self.NONEQ:
    x=[sum([sym[m1][m2]*nq[m2] for m2 in range(3)]) for m1 in range(3)]
    for h1 in self.pm:
     for k1 in self.pm:
      for l1 in self.pm:
       kp=[round(kk,PRECIS) for kk in 
                 (x-(h1*self.e[0]+k1*self.e[1]+l1*self.e[2]))]
       if kp[0]==nq[0] and  kp[1]==nq[1] and  kp[2]==nq[2]:
        found=1
        break
      if found==1: break
     if found==1: break
    if found==1: break
   if found==0: SYMM2.append(sym)

  print(len(SYMM2))
  self.SYMM=SYMM2
 '''
 def calc_weight_of_k(self):
  WK=[0 for i in range(len(self.NONEQ))]
  for i in range(len(self.NONEQ)):
   for j in self.allk:
    if j[3]==i:
     WK[i]+=1
  return WK

 def make_kgrid(self):
  print(' make whole kgrid')
  allk2=[]
  einv=np.linalg.inv(np.transpose(self.e))
  for nq in self.NONEQ:
  # print nq
   for sym in self.SYMM:
    x=[sum([sym[m1][m2]*nq[m2] for m2 in range(3)]) for m1 in range(3)]
    for h1 in self.pm:
     for k1 in self.pm:
      for l1 in self.pm:
       k_point2=[round(kk,PRECIS) for kk in 
                 (x-(h1*self.e[0]+k1*self.e[1]+l1*self.e[2]))]
       k_point3=[round(self.no_of_kpoints[m2]*\
                round(sum([k_point2[m]*einv[m2][m] for m in range(3)]),PRECIS)\
                     ) for m2 in range(3)] #transform from cartesian to crystal coordinates, then we have a cube of points
       if k_point3[0]>=0 and k_point3[0]<=self.no_of_kpoints[0] and \
         k_point3[1]>=0 and k_point3[1]<=self.no_of_kpoints[1] and \
         k_point3[2]>=0 and k_point3[2]<=self.no_of_kpoints[2]:
         allk2.append(k_point2)
         allk2[-1].append(nq[3])
  print (len(allk2))
  allk2=self.sorting(allk2)
  self.allk=self.remove_repeated_items(allk2)

  #rearrange
  allk_in_crystal_coordinates=[ [ round(self.no_of_kpoints[0]*round(sum([v[m]*einv[m2][m] for m in range(3)]),PRECIS)) for m2 in range(3)]+[nv] for nv,v in enumerate(self.allk)]
  allk_in_crystal_coordinates=self.sorting(allk_in_crystal_coordinates)
  self.allk=[ self.allk[i[3]] for i in allk_in_crystal_coordinates if i[0]<self.no_of_kpoints[0]  and i[1]<self.no_of_kpoints[1] and i[2]<self.no_of_kpoints[2]]

  #calc weights of k
  self.WK=self.calc_weight_of_k()

  print(len(self.allk))
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

 def find_k_plus_q(self,k,allk,q):
  einv=np.linalg.inv(np.transpose(self.e))
  kpq=[k[i]+q[i] for i in range(3)]
  found=0
  for sym in self.SYMM:
     x=[sum([sym[m1][m2]*kpq[m2] for m2 in range(3)]) for m1 in range(3)]
     if found==1: break
     for h1 in self.pm:
      if found==1: break
      for k1 in self.pm:
       if found==1: break
       for l1 in self.pm:
        if found==1: break
        k_point2=[round(kk,PRECIS) for kk in 
                  (x-(h1*self.e[0]+k1*self.e[1]+l1*self.e[2]))]
        k_point3=[round(self.no_of_kpoints[m2]*\
                 round(sum([k_point2[m]*einv[m2][m] for m in range(3)]),PRECIS)\
                     ) for m2 in range(3)] #transform from cartesian to crystal coordinates, then we have a cube of points
        if k_point3[0]>=0 and k_point3[0]<=self.no_of_kpoints[0] and \
          k_point3[1]>=0 and k_point3[1]<=self.no_of_kpoints[1] and \
          k_point3[2]>=0 and k_point3[2]<=self.no_of_kpoints[2]:
          kpq2=k_point2
          for ki in allk:
           if found==1: break
           if abs(round(kpq2[0]-ki[0],PRECIS))==0 and abs(round(kpq2[1]-ki[1],PRECIS))==0 and abs(round(kpq2[2]-ki[2],PRECIS))==0:
            kpq_no=ki[3]
            found=1
            break


#  for i in allk:
#   if abs(round(kpq[0]-i[0],PRECIS-3))==0 and abs(round(kpq[1]-i[1],PRECIS-3))==0 and abs(round(kpq[2]-i[2],PRECIS-3))==0:
#    kpq_no=i[3]
#    break
  return [kpq,kpq_no]

