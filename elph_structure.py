import xml.etree.ElementTree as ET
import numpy as np
import structure
PRECIS=6
class elph_structure():
 def __init__(self,ph_structure):
  self.ELPH_sum=[]
  self.elph_dir=ph_structure.elph_dir
  self.prefix=ph_structure.prefix
  self.KPOINTS=[]
  self.WEIGHTS_OF_K=[]
  self.nbnd_el=0
  self.nkp=0
  self.KPOINTS_all=[]
  self.KPOINTS_all_all_q=[]
  self.ALL_COLORS=[] 
  self.lambda_or_elph='' #'lambda' or 'elph'

 def make_kpoints_single_q(self,file,basic_structure):
  self.pm=[0,-1,1]
  tree = ET.parse(self.elph_dir+'elph.'+str(file)+'.1.xml')
  root = tree.getroot()
  self.nkp=int(root.find('PARTIAL_EL_PHON/NUMBER_OF_K').text)
  self.KPOINTS=[ [ round(float(m),PRECIS)\
         for m in root.find(
          'PARTIAL_EL_PHON/K_POINT.'+str(k)+'/COORDINATES_XK'
          ).text.split() ]\
           for k in range(1,self.nkp+1)]
  for numki,ki in enumerate(self.KPOINTS): 
   self.KPOINTS[numki].append(numki)
   found=0
   for h1 in self.pm:
    for k1 in self.pm:
     for l1 in self.pm:
      ki2=[round(kk,PRECIS) for kk in 
                 (ki[0:3]-(h1*basic_structure.e[0]+k1*basic_structure.e[1]+l1*basic_structure.e[2]))]
      for j in basic_structure.allk:
       if ki2[0]==j[0] and ki2[1]==j[1] and ki2[2]==j[2]:
        self.KPOINTS[numki].append(j[3])
        found=1
        break
      if found==1: break
     if found==1: break
    if found==1: break
  structure_new=structure.structure()
  structure_new.NONEQ=self.KPOINTS
  structure_new.SYMM=basic_structure.SYMM
  structure_new.e=basic_structure.e
  structure_new.no_of_kpoints=basic_structure.no_of_kpoints
 # structure_new.check_symm()
  structure_new.make_kgrid()
  self.KPOINTS_all=structure_new.allk
  self.WEIGHTS_OF_K=structure_new.WK
  self.KPOINTS_all_all_q.append(structure_new.allk)


 def w0gauss(self,x):
  degauss=0.02
  sqrtpm1= 1. / 1.77245385090551602729
  # cold smearing  (Marzari-Vanderbilt-DeVita-Payne)
  arg = min (200., (x - 1.0 /  (2.0)**0.5 ) **2.)
  return sqrtpm1 * np.exp ( - arg) * (2.0 - ( 2.0)**0.5 * x)/0.02

 

 def elph_matrix_to_gep(self,nat,dyn,el_ph_mat,w2,pat):
  gep=[[ [[complex(0,0) for ii in range(3*nat)] for ik in range(self.nkp)] for ib in range(self.nbnd_el)] for jb in range(self.nbnd_el)]
  gep2=[[ [[[complex(0,0) for jj in range(3*nat)] for ii in range(3*nat)] for ik in range(self.nkp)] for ib in range(self.nbnd_el)] for jb in range(self.nbnd_el)]
  '''
  for ik in range(self.nkp):
     for ib in range(self.nbnd_el):
      for jb in range(self.nbnd_el):
        for ii in range(3*nat):
          gep[jb][ib][ik][ii] = np.dot(pat[ii], el_ph_mat[jb][ib][ik])
        gep[jb][ib][ik] = np.matmul(gep[jb][ib][ik], dyn)
  '''
  for ik in range(self.nkp):
     for ib in range(self.nbnd_el):
      for jb in range(self.nbnd_el):
        for ii in range(3*nat):
         for jj in range(3*nat):
          gep2[jb][ib][ik][ii][jj] = np.conjugate(el_ph_mat[jb][ib][ik][ii])*el_ph_mat[jb][ib][ik][jj]
        for nu in range(3*nat):
         for mu in range(3*nat):
          for vu in range(3*nat):
           gep[jb][ib][ik][nu] += np.conjugate(dyn[mu][nu])*gep2[jb][ib][ik][mu][vu]*dyn[vu][nu]

  for ik in range(self.nkp):
   for ii in range(3*nat):
    for ib in range(self.nbnd_el):
     for jb in range(self.nbnd_el):
      if w2[ii]<=0.:
        gep[jb][ib][ik][ii]=0.
      else:
        gep[jb][ib][ik][ii]=el_ph_mat[jb][ib][ik][ii]/((w2[ii]**0.5) * 2.)*0.5
  return gep



 def read_elph_single_q(self,file,ph_structure,el_structure): #read info from first file
  tree = ET.parse(self.elph_dir+'elph.'+str(file)+'.1.xml')
  root = tree.getroot()
  self.nbnd_el=int(root.find('PARTIAL_EL_PHON/NUMBER_OF_BANDS').text)
  ELPH=[[[ [] for k in range(self.nkp)] for j in range(self.nbnd_el)] for i in range(self.nbnd_el)] #stores elph[k][ibnd][jbnd][nmode]
  self.ELPH_sum=[[ [] for k in range(self.nkp)] for j in range(self.nbnd_el) ]
#stores elph[k][ibnd][nmode][mmode]

  #read elph  from all files
  for mode in range(1,len(ph_structure.NONDEG[file-1])+1):
   #elph
   tree = ET.parse(self.elph_dir+'elph.'+str(file)+'.'+str(mode)+'.xml')
   root = tree.getroot()
   for country in root.iter('PARTIAL_EL_PHON'):
    for k in range(1,self.nkp+1):
     for town in country.iter('K_POINT.'+str(k)):
      partial_elph=town.find('PARTIAL_ELPH')
      if k==1: npert=int(partial_elph.get('size'))/self.nbnd_el/self.nbnd_el
      elph_k=[ complex(float(m.replace(',',' ').split()[0]), 
                float(m.replace(',',' ').split()[1])) 
                for m in partial_elph.text.split('\n') if len(m.split())>0  ]
      for jband in range(self.nbnd_el):
       for iband in range(self.nbnd_el):
        for iipert in range(npert):
         ELPH[jband][iband][k-1].append\
          (elph_k[jband*self.nbnd_el*npert+iband*npert+iipert])
     
  ELPH=self.elph_matrix_to_gep(ph_structure.nat,\
       ph_structure.DYN[file-1],ELPH,ph_structure.FREQ[file-1],\
       ph_structure.PATT[file-1])

  for jband in range(self.nbnd_el):
    for k in range(1,self.nkp+1):
       self.ELPH_sum[jband][k-1]=[complex(0,0) for i in range(ph_structure.no_of_modes)]
       for iipert in range(ph_structure.no_of_modes):
         for iband in range(self.nbnd_el):
           self.ELPH_sum[jband][k-1][iipert]+=\
                self.w0gauss(el_structure.ef-el_structure.ENE[jband][self.KPOINTS[k-1][4]])*ELPH[jband][iband][k-1][iipert]*self.WEIGHTS_OF_K[k-1]

 def elph_single_q_in_whole_kgrid(self,q,structure,ph_structure,el_structure,l_or_gep):
  self.lambda_or_elph=l_or_gep
 #choose only bands which cross EF and sum over j in pairs <i,j>
  print('From all '+str(self.nbnd_el)+' bands detected in elph calc. only bands ',el_structure.bands_num,' cross EF and will be written in frmsf')
  COLORS=[ [] for jband in el_structure.bands_num]
  if self.lambda_or_elph=='elph':
   for numk,k in enumerate(self.KPOINTS):
    for numjband,jband in enumerate(el_structure.bands_num): 
     COLORS[numjband].append( sum([elph.real**2+elph.imag**2 for elph in self.ELPH_sum[jband][numk]])   )
  else:
   for numk,k in enumerate(self.KPOINTS):
    for numjband,jband in enumerate(el_structure.bands_num): 
     COLORS[numjband].append( sum([(elph.real**2+elph.imag**2)/ph_structure.FREQ[q-1][num] for num,elph in enumerate(self.ELPH_sum[jband][numk])])   )
  
  print len(COLORS),[len(i) for i in COLORS]
  print len(el_structure.ENE_fs),[len(i) for i in el_structure.ENE_fs]
  print len(structure.allk)
  if self.lambda_or_elph=='elph':
   h=open('elph'+str(q)+'.frmsf','w')
  else:
   h=open('lambda'+str(q)+'.frmsf','w')
  h.write(str(structure.no_of_kpoints[0])+' '+str(structure.no_of_kpoints[1])+' ' +str(structure.no_of_kpoints[2])+'\n')
  h.write('1\n'+str(len(el_structure.bands_num))+'\n')
  for i in structure.e:
   for j in i:
    h.write(str(j)+' ')
   h.write('\n')
  for bnd in el_structure.ENE_fs:
   for k in structure.allk:
    h.write(str(bnd[k[3]])+'\n')
  for bnd in COLORS:
   for k in self.KPOINTS_all:
    h.write(str(bnd[k[3]])+'\n')
  h.close()
  self.ALL_COLORS.append(COLORS)


 def sum_over_q(self,ph_structure,structure,el_structure):
  print('summing over q...')
  SUMMED_COLORS=[[0 for k in range(len(self.KPOINTS_all))] for jband in self.ALL_COLORS[0]]
  for band in range(len(SUMMED_COLORS)):
   for numq,col_q in enumerate(self.ALL_COLORS):
    for numk,k in enumerate(self.KPOINTS_all_all_q[numq]):
     SUMMED_COLORS[band][numk]+=col_q[band][k[3]]*ph_structure.multiplicity_of_qs[numq]
  if self.lambda_or_elph=='elph':
   h=open('elph.frmsf','w')
  else:
   h=open('lambda.frmsf','w')
  h.write(str(structure.no_of_kpoints[0])+' '+str(structure.no_of_kpoints[1])+' ' +str(structure.no_of_kpoints[2])+'\n')
  h.write('1\n'+str(len(el_structure.bands_num))+'\n')
  for i in structure.e:
   for j in i:
    h.write(str(j)+' ')
   h.write('\n')
  for bnd in el_structure.ENE_fs:
   for k in structure.allk:
    h.write(str(bnd[k[3]])+'\n')
  for bnd in SUMMED_COLORS:
   for k in bnd:
    h.write(str(k)+'\n')
  h.close()

