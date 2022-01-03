import xml.etree.ElementTree as ET
import numpy as np
import structure,ph_structure
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interpn
from multiprocessing import Process,Pool
PRECIS=6
def symmetrize(a):
    """
    Return a symmetrized version of NumPy array a.

    Values 0 are replaced by the array value at the symmetric
    position (with respect to the diagonal), i.e. if a_ij = 0,
    then the returned array a' is such that a'_ij = a_ji.

    Diagonal values are left untouched.

    a -- square NumPy array, such that a_ij = 0 or a_ji = 0, 
    for i != j.
    """
    return a + a.T - np.diag(a.diagonal())

class elph_structure():
 def __init__(self,phh_structure,lambda_or_elph):
  self.ALL_COLORS=[[] for i in range(len(phh_structure.Q))] 
  self.KPOINTS_all_all_q=[]
  self.SUMMED_COLORS=[]
  self.lambda_or_elph=lambda_or_elph
 def parallel_job(self,structure,el_structure,phh_structure):
  nq=len(phh_structure.Q)
  no_of_pool=nq
  with Pool(no_of_pool) as pol:
   results=pol.map(self.single_job,
              [[int(q+1),structure,phh_structure,el_structure,self.lambda_or_elph] for q in range(nq)])
   self.ALL_COLORS=[ i[0] for i in results]
   self.KPOINTS_all_all_q=[ i[1] for i in results]
 def single_job(self,args):
  [q,structure,phh_structure,el_structure,lambda_or_elph]=args
  print('calculations for '+str(q)+'. of total '+str(len(phh_structure.Q))+' q points')
  elph_q=elph_structure_single_q(phh_structure)
  elph_q.make_kpoints_single_q(q,structure,phh_structure.Q[q-1])
  elph_q.read_elph_single_q(q,phh_structure,el_structure,structure )
  elph_q.elph_single_q_in_whole_kgrid(q,structure,\
                phh_structure,el_structure,lambda_or_elph)  #'lambda' or 'elph'
  return [elph_q.COLORS,elph_q.KPOINTS_all]
  print(str(q)+'. point ended')
  
 def extrapolate_values(self,structure,el_structure):
#  for nk,k in enumerate(structure.allk):
#   str.
  [nx,ny,nz]=structure.no_of_kpoints
  no_of_bands=len(el_structure.ENE_fs)
  ENE=np.zeros(shape=(no_of_bands,nx+1,ny+1,nz+1))
  COL=np.zeros(shape=(no_of_bands,nx+1,ny+1,nz+1))
  for nbnd in range(no_of_bands):
   for x in range(nx+1):
    for y in range(ny+1):
     for z in range(nz+1):
      ENE[nbnd][x,y,z]=el_structure.ENE_fs[nbnd][structure.allk[(x%nx)*ny*ny+(y%ny)*nz+(z%nz)][3]]
      COL[nbnd][x,y,z]=self.SUMMED_COLORS[nbnd][(x%nx)*ny*ny+(y%ny)*nz+(z%nz)]
  '''
  EXTRA_ENE=[ [[[0 for z in range(2*nz)] for y in range(2*ny)] for x in range(2*nx)]  for i in range(len(el_structure.ENE_fs))]
  EXTRA_COLOR=[ [[[0 for z in range(2*nz)] for y in range(2*ny)] for x in range(2*nx)]  for i in range(len(el_structure.ENE_fs))]
  for nbnd,bnd in enumerate(ENE):
   for x in range(2*nx):
    for y in range(2*ny):
     for z in range(2*nz):
      EXTRA_ENE[nbnd][x][y][z]=interpn(mgrid, ENE[nbnd], np.array([x,y,z]))[0]
      EXTRA_COLOR[nbnd][x][y][z]=interpn(mgrid, COL[nbnd], np.array([x,y,z]))[0]
#      self.extrapolate_value(x,y,z,nx,ny,nz,EXTRA_ENE[nbnd],bnd,structure.allk)
#      self.extrapolate_value(x,y,z,nx,ny,nz,EXTRA_COLOR[nbnd],self.SUMMED_COLORS[nbnd],structure.allk)
  '''
  self.X,self.Y,self.Z=np.linspace(0,nx,nx+1),np.linspace(0,ny,ny+1),np.linspace(0,nz,nz+1)
  X2,Y2,Z2=np.linspace(0,nx-.5,2*nx),np.linspace(0,ny-.5,2*ny),np.linspace(0,nz-.5,2*nz)
  ene_interp=[RegularGridInterpolator((self.X,self.Y,self.Z), bnd) for bnd in ENE]
  col_interp=[RegularGridInterpolator((self.X,self.Y,self.Z), bnd) for bnd in ENE]
  EXTRA_ENE=[[[[ene_interp[nbnd]([i,j,k])[0]  for k in Z2 ] for j in Y2 ] for i in X2 ] for nbnd in range(len(ENE))]
  EXTRA_COL=[[[[col_interp[nbnd]([i,j,k])[0]  for k in Z2 ] for j in Y2 ] for i in X2 ] for nbnd in range(len(COL))]
  h=open('lambda_extrapolated.frmsf','w')
  h.write(str(2*nx)+' '+str(2*ny)+' ' +str(2*nz)+'\n')
  h.write('1\n'+str(len(el_structure.bands_num))+'\n')
  for i in structure.e:
   for j in i:
    h.write(str(j)+' ')
   h.write('\n')
  for bnd in EXTRA_ENE:
   for x in range(2*nx):
    for y in range(2*ny):
     for z in range(2*nz):
      h.write(str(bnd[x][y][z])+'\n')
  for bnd in EXTRA_COL:
   for x in range(2*nx):
    for y in range(2*ny):
     for z in range(2*nz):
      h.write(str(bnd[x][y][z])+'\n')

 def sum_over_q(self,phh_structure,structure,el_structure):
  print('summing over q...')
  self.SUMMED_COLORS=[[0 for k in range(len(structure.allk))] for jband in range(el_structure.fermi_nbnd_el)]
  for band in range(len(self.SUMMED_COLORS)):
   for numq,col_q in enumerate(self.ALL_COLORS):
    for numk,k in enumerate(structure.allk):
     self.SUMMED_COLORS[band][numk]+=col_q[band][numk] #*phh_structure.multiplicity_of_qs[numq]
  print ('summed lambda=',sum([ sum(i) for i in self.SUMMED_COLORS]))
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
  for bnd in self.SUMMED_COLORS:
   for k in bnd:
    h.write(str(k)+'\n')
  h.close()

 '''
 def extrapolate_value(self,x,y,z,nx,ny,nz,bnd,allk):
  suma=0
  nsuma=0
  if x!=0: 
   kp1=x*ny*ny+y*nz+z
   suma+=bnd[allk[kp1][3]]
   nsuma+=1
  if x!=nx-1: 
   kp1=(x+1)*ny*ny+y*nz+z
   suma+=bnd[allk[kp1][3]]
   nsuma+=1
  if y!=0: 
   kp1=x*ny*ny+y*nz+z
   suma+=bnd[allk[kp1][3]]
   nsuma+=1
  if y!=ny-1: 
   kp1=x*ny*ny+(y+1)*nz+z
   suma+=bnd[allk[kp1][3]]
   nsuma+=1
  if z!=0: 
   kp1=x*ny*ny+y*nz+z-1
   suma+=bnd[allk[kp1][3]]
   nsuma+=1
  if z!=nz-1: 
   kp1=x*ny*ny+y*nz+z+1
   suma+=bnd[allk[kp1][3]]
   nsuma+=1
  if suma!=0: return suma/nsuma
  else: return 0
  '''



  
class elph_structure_single_q():
 def __init__(self,phh_structure):
  self.ELPH_sum=[]
  self.elph_dir=phh_structure.elph_dir
  self.prefix=phh_structure.prefix
  self.KPOINTS=[] #[0:2] list of noneq k at given q [3] its no in list of nonequivalent kpoints at given q  (self.KPOINTS)  and [4] its no in list of nonequivalent kpoints of electronic structure (el_structure.NONEQ)
  self.WEIGHTS_OF_K=[]
  self.nbnd_el=0
  self.fermi_nbnd_el=0
  self.nkp=0
  self.KPOINTS_all=[] #[0:2] list of allk [3] its no in list of nonequivalent kpoints at given q  (self.KPOINTS)  
  self.KPQ=[] #[0] k+q, [1] its no in list of nonequivalent kpoints at given q  (self.KPOINTS)  
  self.lambda_or_elph='' #'lambda' or 'elph'
  self.COLORS=[]

 def make_kpoints_single_q(self,q_no,basic_structure,q):
  print('make kgrid at given q')
  self.pm=[0,-1,1]
  tree = ET.parse(self.elph_dir+'elph.'+str(q_no)+'.1.xml')
  root = tree.getroot()
  self.nkp=int(root.find('PARTIAL_EL_PHON/NUMBER_OF_K').text)
  self.KPOINTS=[ [ round(float(m),PRECIS)\
         for m in root.find(
          'PARTIAL_EL_PHON/K_POINT.'+str(k)+'/COORDINATES_XK'
          ).text.split() ]\
           for k in range(1,self.nkp+1)]
  for numki,ki in enumerate(self.KPOINTS): 
   self.KPOINTS[numki].append(numki)
   '''
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
    '''
  structure_new=structure.structure()
  structure_new.NONEQ=list(self.KPOINTS)
  structure_new.SYMM=list(basic_structure.SYMM)
  structure_new.SYMM_crystal=list(basic_structure.SYMM_crystal)
  structure_new.e=basic_structure.e
  structure_new.no_of_kpoints=basic_structure.no_of_kpoints
  structure_new.check_symm(q,basic_structure.NONEQ)
  print(q_no,'no of sym',len(structure_new.SYMM))
  self.SYMM_q=structure_new.SYMM_crystal
  structure_new.make_kgrid(q_no)
  structure_new.find_newkpoints_in_old_list(basic_structure.allk) #aatach 4-th element to each k:  #it's nonequivalent in basic k-point list. Added to both allk AND NONEQ
  self.KPOINTS_all=structure_new.allk
  self.WEIGHTS_OF_K=structure_new.WK
  self.KPOINTS=structure_new.NONEQ
#  print(q_no,[ i for i in self.KPOINTS if [4]==None])
#  self.KPOINTS_all_all_q.append(structure_new.allk)

  self.KPQ=[]
  for k in self.KPOINTS:
   self.KPQ.append(structure_new.find_k_plus_q(k, self.KPOINTS_all,q)) #[kpq,no its nonequivalent in kpoints list for this q, no its nonequivalent in basic kpoint list]


 def w0gauss(self,x):
  degauss=0.02
  x2=x/degauss
  sqrtpm1= 1. / 1.77245385090551602729
  # cold smearing  (Marzari-Vanderbilt-DeVita-Payne)
  arg = min (200., (x2 - 1.0 /  (2.0)**0.5 ) **2.)
  return sqrtpm1 * np.exp ( - arg) * (2.0 - ( 2.0)**0.5 * x2)/degauss

 
 def read_elph_single_q(self,q_point_no,phh_structure,el_structure,structure): 
  print(str(q_point_no)+': read elph from file...')
  #read info from first file
  tree = ET.parse(self.elph_dir+'elph.'+str(q_point_no)+'.1.xml')
  root = tree.getroot()
  self.nbnd_el=int(root.find('PARTIAL_EL_PHON/NUMBER_OF_BANDS').text)
  #choose only bands which cross EF and sum over j in pairs <i,j>
  print(str(q_point_no)+': From all '+str(self.nbnd_el)+' bands detected in elph calc. only bands ',el_structure.bands_num,' cross EF and will be written in frmsf')
  self.fermi_nbnd_el=len(el_structure.bands_num)
#  ELPH=[[[ [] for k in range(self.nkp)] for j in range(self.nbnd_el)] for i in range(self.fermi_nbnd_el)] #stores elph[k][ibnd][jbnd][nmode]
  ELPH=np.zeros(shape=(self.fermi_nbnd_el,self.fermi_nbnd_el,\
                self.nkp,phh_structure.no_of_modes), dtype=complex) #stores elph[k][ibnd][jbnd][nmode]
  ELPH2=np.zeros(shape=(phh_structure.no_of_modes,phh_structure.no_of_modes), dtype=complex) #stores elph[k][ibnd][jbnd][nmode]
#  self.ELPH_sum1=np.zeros(shape=(self.fermi_nbnd_el,\
#           self.nkp,phh_structure.no_of_modes,phh_structure.no_of_modes),dtype=complex)
  self.ELPH_sum=np.zeros(shape=(self.fermi_nbnd_el,\
           self.nkp,phh_structure.no_of_modes),dtype=complex)
#  self.ELPH_sum0=np.zeros(shape=(self.fermi_nbnd_el,\
#           self.nkp,phh_structure.no_of_modes,phh_structure.no_of_modes),dtype=complex)
#stores elph[ibnd][k][nmode][mmode]

  #read elph  from all files
  imode=0
  for mode in range(1,len(phh_structure.NONDEG[q_point_no-1])+1):
   #elph
   tree = ET.parse(self.elph_dir+'elph.'+str(q_point_no)+'.'+str(mode)+'.xml')
   root = tree.getroot()
   for country in root.iter('PARTIAL_EL_PHON'):
    if int(country.find('NUMBER_OF_BANDS').text)!=self.nbnd_el: print(str(q_point_no)+'Warning. No of bands!= no of bands from el structure')
    for k in range(1,self.nkp+1):
     for town in country.iter('K_POINT.'+str(k)):
      partial_elph=town.find('PARTIAL_ELPH')
      if k==1: 
       npert0=int(partial_elph.get('size'))/self.nbnd_el/self.nbnd_el
       npert=int(npert0)
       print( npert, mode,phh_structure.ORDER_OF_IRR[q_point_no-1],phh_structure.FREQ[q_point_no-1])
       if npert/npert0!=1: 
        print(str(q_point_no)+"WARNING: npert is not int, but is equal to"+str(npert0))
      elph_k=[ complex(float(m.replace(',',' ').split()[0]), 
                       float(m.replace(',',' ').split()[1])) 
                for m in partial_elph.text.split('\n') if len(m.split())>0  ]
      for iipert in range(npert):
       nmode=imode+iipert #phh_structure.ORDER_OF_IRR[q_point_no-1][imode+iipert]
#       elph_k_ii=symmetrize((elph_k_ii))
       for numjband,jband in enumerate(el_structure.bands_num):
        for numiband,iband in enumerate(el_structure.bands_num):
          ELPH[numjband][numiband][k-1][nmode]=elph_k[iipert*self.nbnd_el*npert+jband*self.nbnd_el+iband]
#          elph_k[jband*self.nbnd_el*npert+iband*npert+iipert] #--wrong
#          elph_k[iipert*self.nbnd_el*npert+iband*self.nbnd_el+jband] #--rather wrong
#              CALL iotk_write_dat(iunpun, "PARTIAL_ELPH", &
#                                         el_ph_mat_rec_col(:,:,ik,:))
#		el_ph_mat_rec_col(nbnd,nbnd,nksqtot,npe)
   imode=imode+npert   

  
  for iband in range(self.fermi_nbnd_el):
#   sumweight=0
   for k in range(1,self.nkp+1):
    ELPH2=np.zeros(shape=(phh_structure.no_of_modes,phh_structure.no_of_modes), dtype=complex) #stores elph for given iband, k [summed over jband
    for jband in range(self.fermi_nbnd_el):
     weight=self.w0gauss(-el_structure.ENE_fs[jband][self.KPQ[k-1][2]])*self.WEIGHTS_OF_K[self.KPQ[k-1][1]]
#     sumweight+=weight
     for iipert in range(phh_structure.no_of_modes):
      for jjpert in range(phh_structure.no_of_modes):
       ELPH2[iipert][jjpert]+=np.conjugate(ELPH[jband][iband][k-1][iipert])*ELPH[jband][iband][k-1][jjpert]*weight
    ELPH2=ph_structure.symmetrize(phh_structure.nat,np.array(phh_structure.PATT[q_point_no-1]),
ELPH2, structure.at,structure.e, structure.SYMM_crystal, self.SYMM_q)
    for nu in range(phh_structure.no_of_modes):
      for mu in range(phh_structure.no_of_modes):
       for vu in range(phh_structure.no_of_modes):
         self.ELPH_sum[iband][k-1][nu] += np.conjugate(phh_structure.DYN[q_point_no-1][0][nu][mu])*ELPH2[vu][mu]*phh_structure.DYN[q_point_no-1][0][nu][vu]     

    #print (self.ELPH_sum[iband][k-1])

  '''
  #złes
  for k in range(1,self.nkp+1):
   for iband in range(self.fermi_nbnd_el):
    for jband in range(self.fermi_nbnd_el):
     weight=self.w0gauss(-el_structure.ENE_fs[jband][self.KPQ[k-1][2]])*self.WEIGHTS_OF_K[self.KPQ[k-1][1]]
     for iipert in range(phh_structure.no_of_modes):
       self.ELPH_sum[iband][k-1][iipert] += np.conjugate(ELPH[iband][jband][k-1][iipert])*ELPH[iband][jband][k-1][iipert]*weight
  '''
  '''
  #sum over jband
   for iband in range(self.fermi_nbnd_el):
      sum_weight=0
      for iipert in range(phh_structure.no_of_modes):
       w2ii=phh_structure.FREQ[q_point_no-1][iipert]
       for jjpert in range(phh_structure.no_of_modes):
        w2jj=phh_structure.FREQ[q_point_no-1][jjpert]
        if w2ii>=0.0001 or w2jj>=0.0001:
         for jband in range(self.nbnd_el):
           weight=self.w0gauss(-el_structure.ENE[jband][self.KPOINTS[self.KPQ[k-1][1]][4]])*self.WEIGHTS_OF_K[self.KPQ[k-1][1]]
           gep2 = np.conjugate(ELPH[jband][iband][k-1][iipert])\
		  *ELPH[jband][iband][k-1][jjpert] \
                  /(w2ii**0.5)/(w2jj**0.5)
           self.ELPH_sum[iband][k-1]+=weight*gep2
           sum_weight+=weight
      if sum_weight!=0: self.ELPH_sum[iband][k-1]=self.ELPH_sum[iband][k-1]/sum_weight
  '''
  '''
  #sum over jband
  #GORZEJ
  for iband in range(self.fermi_nbnd_el):
    for k in range(1,self.nkp+1):
      sum_weight=0
      for iipert in range(phh_structure.no_of_modes):
       w2ii=phh_structure.FREQ[q_point_no-1][iipert]
       if w2ii>=0.0001:
        for jband in range(self.nbnd_el):
           weight=self.w0gauss(-el_structure.ENE[jband][self.KPOINTS[self.KPQ[k-1][1]][4]])*self.WEIGHTS_OF_K[self.KPQ[k-1][1]]
           gep2 = np.conjugate(ELPH[jband][iband][k-1][iipert])\
		  *ELPH[jband][iband][k-1][iipert]\
                  /w2ii
           self.ELPH_sum[iband][k-1]+=weight*gep2
           sum_weight+=weight
      if sum_weight!=0: self.ELPH_sum[iband][k-1]=self.ELPH_sum[iband][k-1]/sum_weight
  '''

 def elph_single_q_in_whole_kgrid(self,q,structure,phh_structure,el_structure,l_or_gep):
  print(str(q)+': print elph to file')
  self.lambda_or_elph=l_or_gep
  self.COLORS=np.zeros(shape=(self.fermi_nbnd_el, self.nkp))
  print(max([i[4] for i in self.KPOINTS]),[len(bnd) for bnd in el_structure.ENE_fs])
  wagi0=[[self.w0gauss(-el_structure.ENE_fs[jband][k[4]]) for numk,k in enumerate(self.KPOINTS)] for jband in range(len(self.COLORS)) ]
  wagi=[[wagi0[jband][numk]*self.WEIGHTS_OF_K[numk] for numk in range(len(self.KPOINTS))] for jband in range(len(self.COLORS)) ]
  waga_sum=sum([sum(w) for w in wagi])
  if self.lambda_or_elph=='elph':
   for numk,k in enumerate(self.KPOINTS):
    for jband in range(self.fermi_nbnd_el): 
     self.COLORS[jband][numk]= np.abs(np.sum(self.ELPH_sum[jband][numk]))  
  else:
   for numk,k in enumerate(self.KPOINTS):
    for jband in range(self.fermi_nbnd_el): 
     self.COLORS[jband][numk]= 2.*np.abs(np.sum([elph/(phh_structure.FREQ[q-1][phh_structure.ORDER_OF_IRR[q-1][num]]) for num,elph in enumerate(self.ELPH_sum[jband][numk])])) #*wagi0[jband][numk]/waga_sum  #abs=sqrt(real^2+im^2)
#  print (q,': lambda_sum=',lambda_sum)
  lambda_sum=sum([sum([ self.COLORS[jband][numk]*self.WEIGHTS_OF_K[numk] for numk in range(len(self.COLORS[jband]))]) for jband in range(len(self.COLORS)) ])
  print (q,'\n\n\n: lambda_sum=',lambda_sum)

  self.COLORS=[[ bnd[k[3]] for k in self.KPOINTS_all] for bnd in self.COLORS]


  print (q,":",len(self.COLORS),[len(i) for i in self.COLORS])
  print (q,":",len(el_structure.ENE_fs),[len(i) for i in el_structure.ENE_fs])
  print (q,":",len(structure.allk))
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
  for bnd in self.COLORS:
   for numk,k in enumerate(self.KPOINTS_all):
    h.write(str(bnd[numk])+'\n')
  h.close()
#  self.ALL_COLORS[q-1]=COLORS





 def extrapolate_value(self,x,y,z,nx,ny,nz,ENE,bnd,allk):
      kp=x*ny*ny+  y*nz  +   z
      ENE[2*x][2*y][2*z]=bnd[allk[kp][3]]
      if x!=nx-1: 
       divby=2
       ENE[2*x+1][2*y][2*z]=(bnd[allk[kp][3]]+bnd[allk[(x+1)*ny*ny+y*nz+z][3]]) #/2.
       if y!=0: 
        ENE[2*x+1][2*y][2*z]+=ENE[2*x+1][2*y-1][2*z]
        divby+=1
       if z!=0: 
        ENE[2*x+1][2*y][2*z]+=ENE[2*x+1][2*y][2*z-1]
        divby+=1
       ENE[2*x+1][2*y][2*z]=ENE[2*x+1][2*y][2*z]/divby

      if y!=ny-1: 
       divby=2
       ENE[2*x][2*y+1][2*z]=(bnd[allk[kp][3]]+bnd[allk[x*ny*ny+(y+1)*nz+z][3]]) #/2.
       if x!=0: 
        ENE[2*x][2*y+1][2*z]+=ENE[2*x-1][2*y+1][2*z]
        divby+=1
       if z!=0: 
        ENE[2*x][2*y+1][2*z]+=ENE[2*x][2*y+1][2*z-1]
        divby+=1
       ENE[2*x][2*y+1][2*z]=ENE[2*x][2*y+1][2*z]/divby

      if z!=nz-1: 
       divby=2
       ENE[2*x][2*y][2*z+1]=(bnd[allk[kp][3]]+bnd[allk[x*ny*ny+y*nz+z+1][3]]) #/2.
       if x!=0: 
        ENE[2*x][2*y][2*z+1]+=ENE[2*x-1][2*y][2*z+1]
        divby+=1
       if y!=0: 
        ENE[2*x][2*y][2*z+1]+=ENE[2*x][2*y-1][2*z+1]
        divby+=1
       ENE[2*x][2*y][2*z+1]=ENE[2*x][2*y][2*z+1]/divby


      if x!=nx-1 and y!=ny-1: 
       ENE[2*x+1][2*y+1][2*z]=\
        (bnd[allk[kp][3]]+bnd[allk[(x+1)*ny*ny+(y+1)*nz+z][3]]\
        +bnd[allk[(x+1)*ny*ny+y*nz+z][3]]+bnd[allk[x*ny*ny+(y+1)*nz+z][3]]\
        )/4.
       if z!=0: ENE[2*x+1][2*y+1][2*z]=(ENE[2*x+1][2*y+1][2*z]*4+ENE[2*x+1][2*y+1][2*z-1]+ENE[2*x+1][2*y][2*z-1]+ENE[2*x][2*y+1][2*z-1])/7

      if y!=ny-1 and z!=nz-1: 
       ENE[2*x][2*y+1][2*z+1]=\
        (bnd[allk[kp][3]]+bnd[allk[x*ny*ny+(y+1)*nz+z+1][3]]\
        +bnd[allk[x*ny*ny+(y+1)*nz+z][3]]+bnd[allk[x*ny*ny+y*nz+z+1][3]]\
        )/4.
       if x!=0: ENE[2*x][2*y+1][2*z+1]=(ENE[2*x][2*y+1][2*z+1]*4+ENE[2*x-1][2*y+1][2*z+1])/5

      if x!=nx-1 and z!=nz-1: 
       ENE[2*x+1][2*y][2*z+1]=\
        (bnd[allk[kp][3]]+bnd[allk[(x+1)*ny*ny+y*nz+z+1][3]]\
        +bnd[allk[(x+1)*ny*ny+y*nz+z][3]]+bnd[allk[x*ny*ny+y*nz+z+1][3]]\
        )/4.
       if y!=0: ENE[2*x+1][2*y][2*z+1]=(ENE[2*x+1][2*y][2*z+1]*4+ENE[2*x+1][2*y-1][2*z+1])/5

      if x!=nx-1 and y!=ny-1 and z!=nz-1: 
       ENE[2*x+1][2*y+1][2*z+1]=\
        (bnd[allk[kp][3]]+bnd[allk[(x+1)*ny*ny+(y+1)*nz+z+1][3]]\
        +bnd[allk[(x+1)*ny*ny+y*nz+z+1][3]]+bnd[allk[x*ny*ny+(y+1)*nz+z+1][3]] +bnd[allk[(x+1)*ny*ny+(y+1)*nz+z][3]]\
        +bnd[allk[(x+1)*ny*ny+y*nz+z][3]]+bnd[allk[x*ny*ny+(y+1)*nz+z][3]] +bnd[allk[x*ny*ny+y*nz+z+1][3]]\
        )/8.



