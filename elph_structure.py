import xml.etree.ElementTree as ET
import numpy as np
import structure,ph_structure,el_structure
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interpn
from multiprocessing import Process,Pool
import copy
PRECIS=6
RY_TO_THZ=3289.8449

#def symmetrize(a):
"""
    Return a symmetrized version of NumPy array a.

    Values 0 are replaced by the array value at the symmetric
    position (with respect to the diagonal), i.e. if a_ij = 0,
    then the returned array a' is such that a'_ij = a_ji.

    Diagonal values are left untouched.

    a -- square NumPy array, such that a_ij = 0 or a_ji = 0, 
    for i != j.
"""
 #   return a + a.T - np.diag(a.diagonal())
#
class elph_structure():
 def __init__(self,phh_structure,lambda_or_elph):
  self.ALL_COLORS=[[] for i in range(len(phh_structure.Q))] 
  self.KPOINTS_all_all_q=[]
  self.SUMMED_COLORS=[]
  self.lambda_or_elph=lambda_or_elph
 def parallel_job(self,structure,ell_structure,phh_structure):
  nq=len(phh_structure.Q)
  no_of_pool=nq
  with Pool(no_of_pool) as pol:
   results=pol.map(self.single_job,
              [[int(q+1),structure,phh_structure,ell_structure,self.lambda_or_elph] for q in range(nq)])
   self.ALL_COLORS=[ i[0] for i in results]
   self.KPOINTS_all_all_q=[ i[1] for i in results]
 def single_job(self,args):
  [q,structure,phh_structure,ell_structure,lambda_or_elph]=args
  print('calculations for '+str(q)+'. of total '+str(len(phh_structure.Q))+' q points')
  elph_q=elph_structure_single_q(q,phh_structure)
  elph_q.make_kpoints_single_q(q,structure,phh_structure.Q[q-1],phh_structure.Q_crystal[q-1],phh_structure.qstar[q-1])
  elph_q.read_elph_single_q(q,phh_structure,ell_structure,structure )
  elph_q.elph_single_q_in_whole_kgrid(q,structure,\
                phh_structure,ell_structure,lambda_or_elph)  #'lambda' or 'elph'
  return [elph_q.COLORS,elph_q.KPOINTS_all]
  print(str(q)+'. point ended')
  


 def sum_over_q(self,phh_structure,structure,ell_structure):
  print('summing over q...')
  self.SUMMED_COLORS=[[0 for k in range(len(structure.allk))] for jband in range(ell_structure.fermi_nbnd_el)]
#  sum_wag=0
  for band in range(len(self.SUMMED_COLORS)):
   for numk,k in enumerate(structure.allk):
 #   waga=el_structure.w0gauss(-ell_structure.ENE_fs[band][k[3]])
 #   sum_wag+=waga
    for numq,col_q in enumerate(self.ALL_COLORS):
     self.SUMMED_COLORS[band][numk]+=col_q[band][numk] *phh_structure.multiplicity_of_qs[numq] # *waga--done for every q
  print ('summed lambda=',sum([ sum(i) for i in self.SUMMED_COLORS]))
#  print ('summed wagi=',sum_wag)
  if self.lambda_or_elph=='elph':
   h=open('elph.frmsf','w')
  else:
   h=open('lambda.frmsf','w')
  h.write(str(structure.no_of_kpoints[0])+' '+str(structure.no_of_kpoints[1])+' ' +str(structure.no_of_kpoints[2])+'\n')
  h.write('1\n'+str(len(ell_structure.bands_num))+'\n')
  for i in structure.e:
   for j in i:
    h.write(str(j)+' ')
   h.write('\n')
  for bnd in ell_structure.ENE_fs:
   for k in structure.allk:
    h.write(str(bnd[k[3]])+'\n')
  for bnd in self.SUMMED_COLORS:
   for k in bnd:
    h.write(str(np.abs(k))+'\n')
  h.close()



  
class elph_structure_single_q():
 def __init__(self,q,phh_structure):
  self.ELPH_sum=[]
  self.elph_dir=phh_structure.elph_dir
  self.prefix=phh_structure.prefix
  self.no_of_s_q=phh_structure.no_of_s_q[q-1]
  self.KPOINTS=[] #[0:2] list of noneq k at given q [3] its no in list of nonequivalent kpoints at given q  (self.KPOINTS)  and [4] its no in list of nonequivalent kpoints of electronic structure (ell_structure.NONEQ)
  self.WEIGHTS_OF_K=[]
  self.nbnd_el=0
  self.fermi_nbnd_el=0
  self.nkp=0
  self.KPOINTS_all=[] #[0:2] list of allk [3] its no in list of nonequivalent kpoints at given q  (self.KPOINTS)  
  self.KPQ=[] #[0] k+q, [1] its no in list of nonequivalent kpoints at given q  (self.KPOINTS)  
  self.lambda_or_elph='' #'lambda' or 'elph'
  self.COLORS=[]

 def make_kpoints_single_q(self,q_no,basic_structure,q,q_cryst,qstar):
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
  structure_new=copy.deepcopy(basic_structure)
  structure_new.NONEQ=list(self.KPOINTS)
  structure_new.calc_noneq_cryst() #results in NONEQ_cryst
  structure_new.no_of_kpoints=basic_structure.no_of_kpoints
  structure_new.check_symm(q_cryst,basic_structure.NONEQ,self.no_of_s_q,q_no)
  structure_new.calc_irt()
  self.irt=structure_new.irt
  self.rtau=structure_new.rtau
  self.SYMM_q=structure_new.SYMM_crystal
  structure_new.make_kgrid(q_no)
  structure_new.find_newkpoints_in_old_list(basic_structure.allk) #aatach 4-th element to each k:  #it's nonequivalent in basic k-point list. Added to both allk AND NONEQ
  self.KPOINTS_all=structure_new.allk
  self.KPOINTS_all_cryst=structure_new.allk_in_crystal_coordinates
  self.WEIGHTS_OF_K=structure_new.WK
  self.KPOINTS=structure_new.NONEQ
  self.KPOINTS_cryst=structure_new.NONEQ_cryst
#  print(q_no,[ i for i in self.KPOINTS if [4]==None])
#  self.KPOINTS_all_all_q.append(structure_new.allk)

 # print(self.KPOINTS)
  '''
  print('star calculations')
  self.KPQ=[[] for qq in qstar]
  for nqq,qq in enumerate(qstar):
   for k in self.KPOINTS:
    self.KPQ[nqq].append(structure_new.find_k_plus_q(k, self.KPOINTS_all,qq)) #[kpq,no its nonequivalent in kpoints list for this q, no its nonequivalent in basic kpoint list]
  print('star ended')
  '''
  self.KPQ=[]
  for k in  self.KPOINTS_cryst:
   self.KPQ.append(structure_new.find_k_plus_q(k, self.KPOINTS_all,self.KPOINTS_all_cryst,q_cryst))
  h=open('q_'+str(q_no)+'.info','w')
  h.write(str(len(self.SYMM_q))+'\n')
  #for i in self.KPQ:
  # for j in i:
  #  h.write(str(j)+' ')
  # h.write('\n')

  multik=np.zeros((len(basic_structure.NONEQ)))
  for k in self.KPQ: multik[k[2]]+=1
  for ni,i in  enumerate( self.KPQ):
    h.write(str(ni)+' '+str(i)+'\n')
  h.close()
  for i in multik: 
   if i==0: raise ValueError(str(q_no)+' len of multik = 0')

  '''
  #check if kpq are ok
  structure_newq=structure.structure()
  structure_newq.NONEQ=[m[0]+[nm] for nm,m in enumerate(self.KPQ)]
  structure_newq.SYMM=list(structure_new.SYMM)
  structure_newq.SYMM_crystal=list(structure_new.SYMM_crystal)
  structure_newq.e=basic_structure.e
  structure_newq.no_of_kpoints=basic_structure.no_of_kpoints
  structure_newq.make_kgrid(q_no)
  structure_newq.find_newkpoints_in_old_list(basic_structure.allk) #aatach 4-th element to each k:  #it's nonequiv
  self.KPQ_all=structure_newq.allk
  '''

 def read_elph_single_q(self,q_point_no,phh_structure,ell_structure,structure): 
  print(str(q_point_no)+': read elph from file...')
  #read info from first file
  tree = ET.parse(self.elph_dir+'elph.'+str(q_point_no)+'.1.xml')
  root = tree.getroot()
  self.nbnd_el=int(root.find('PARTIAL_EL_PHON/NUMBER_OF_BANDS').text)
  #choose only bands which cross EF and sum over j in pairs <i,j>
  print(str(q_point_no)+': From all '+str(self.nbnd_el)+' bands detected in elph calc. only bands ',ell_structure.bands_num,' cross EF and will be written in frmsf')
  self.fermi_nbnd_el=len(ell_structure.bands_num)
#  ELPH=[[[ [] for k in range(self.nkp)] for j in range(self.nbnd_el)] for i in range(self.fermi_nbnd_el)] #stores elph[k][ibnd][jbnd][nmode]
  ELPH=np.zeros(shape=(self.fermi_nbnd_el,self.fermi_nbnd_el,\
                self.nkp,phh_structure.no_of_modes), dtype=complex) #stores elph[k][ibnd][jbnd][nmode]

  ELPH2=np.zeros(shape=(len(structure.NONEQ),phh_structure.no_of_modes,phh_structure.no_of_modes), dtype=complex) #stores elph[k][ibnd][jbnd][nmode]
#  self.ELPH_sum1=np.zeros(shape=(self.fermi_nbnd_el,\
#           self.nkp,phh_structure.no_of_modes,phh_structure.no_of_modes),dtype=complex)
  self.ELPH_sum=np.zeros(shape=(self.fermi_nbnd_el,\
           len(structure.NONEQ),phh_structure.no_of_modes),dtype=complex)
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
       for numiband,iband in enumerate(ell_structure.bands_num):
        for numjband,jband in enumerate(ell_structure.bands_num):
         ELPH[numiband][numjband][k-1][nmode]=elph_k[iipert*self.nbnd_el*self.nbnd_el+jband*self.nbnd_el+iband]
   #      ELPH[numiband][numjband][k-1][nmode]=elph_k[jband*self.nbnd_el*npert+iband*npert+iipert] #--wrong
   #      ELPH[numiband][numjband][k-1][nmode]=elph_k[iband*self.nbnd_el*npert+iipert*self.nbnd_el+jband] #--wrong
#          elph_k[iipert*self.nbnd_el*npert+iband*self.nbnd_el+jband] #--rather wrong
#              CALL iotk_write_dat(iunpun, "PARTIAL_ELPH", &
#                                         el_ph_mat_rec_col(:,:,ik,:))
#		el_ph_mat_rec_col(nbnd,nbnd,nksqtot,npe)
   imode=imode+npert   



  multik=np.zeros((len(structure.NONEQ)))
  for k in self.KPQ: multik[k[2]]+=1
  for i in multik: 
   if i==0: raise ValueError('len of multik = 0')
  '''
  #checks if kpq are ok
  h=open('FS_'+str(q_point_no)+'.frmsf','w')
  h.write(str(structure.no_of_kpoints[0])+' '+str(structure.no_of_kpoints[1])+' ' +str(structure.no_of_kpoints[2])+'\n')
  h.write('1\n'+str(len(ell_structure.bands_num))+'\n')
  for i in structure.e:
   for j in i:
    h.write(str(j)+' ')
   h.write('\n')
  for bnd in ell_structure.ENE_fs:
   for k in self.KPQ_all:
    h.write(str(bnd[k[4]])+'\n')
  h.close()
  '''

  for iband in range(self.fermi_nbnd_el):
#   sumweight=0
    for numk, k in enumerate(self.KPOINTS):
     weight2=el_structure.w0gauss(-ell_structure.ENE_fs[iband][k[4]])  
     for jband in range(self.fermi_nbnd_el):
      weight=el_structure.w0gauss(-ell_structure.ENE_fs[jband][self.KPQ[numk][2]])  *self.WEIGHTS_OF_K[self.KPQ[numk][1]]
      for iipert in range(phh_structure.no_of_modes):
       for jjpert in range(phh_structure.no_of_modes):
        ELPH2[k[4]][iipert][jjpert]+=np.conjugate(ELPH[iband][jband][numk][iipert])*ELPH[iband][jband][numk][jjpert]*weight*weight2 /multik[k[4]]  
    for k in range(len(structure.NONEQ)):  
     ELPH2[k]=ph_structure.symmetrize(phh_structure.nat,np.array(phh_structure.PATT[q_point_no-1]),
ELPH2[k], structure.at,structure.e, structure.SYMM_crystal, self.SYMM_q,self.irt,self.rtau,phh_structure.Q_crystal[q_point_no-1] )
# dyn=ph_structure.symmetrize(phh_structure.nat,np.array(phh_structure.PATT[qno]),
#     (phh_structure.DYN2[qno]), sstructure.at,sstructure.e,
#     sstructure.SYMM_crystal, structure_new.SYMM_crystal,
#     structure_new.irt,structure_new.rtau,phh_structure.Q_crystal_orig[qno] )
     for nu in range(phh_structure.no_of_modes):
       for mu in range(phh_structure.no_of_modes):
        for vu in range(phh_structure.no_of_modes):
          self.ELPH_sum[iband][k][nu] += np.conjugate(phh_structure.DYN[q_point_no-1][0][nu][mu])*ELPH2[k][vu][mu]*phh_structure.DYN[q_point_no-1][0][nu][vu]  #/multik[k[4]]    

  '''
  multik=np.zeros((len(structure.NONEQ)))
  for k in self.KPOINTS: multik[k[4]]+=1
  for iband in range(self.fermi_nbnd_el):
   for numk, k in enumerate(self.KPOINTS):
      ELPH2=np.zeros(shape=(phh_structure.no_of_modes,phh_structure.no_of_modes), dtype=complex) #stores elph for given iband, k [summed over jband
      for jband in range(self.fermi_nbnd_el):
       weight=el_structure.w0gauss(-ell_structure.ENE_fs[jband][self.KPQ[numk][2]]) *self.WEIGHTS_OF_K[self.KPQ[numk][1]]
       ELPH3=np.zeros(shape=(phh_structure.no_of_modes,phh_structure.no_of_modes), dtype=complex) #stores 
       for iipert in range(phh_structure.no_of_modes):
        for jjpert in range(phh_structure.no_of_modes):
         ELPH3[iipert][jjpert]=np.conjugate(ELPH[iband][jband][numk][iipert])*ELPH[iband][jband][numk][jjpert]
       ELPH3=ph_structure.symmetrize(phh_structure.nat,np.array(phh_structure.PATT[q_point_no-1]),
ELPH3, structure.at,structure.e, structure.SYMM_crystal, self.SYMM_q)
       ELPH2 += ELPH3*weight
      for nu in range(phh_structure.no_of_modes):
        for mu in range(phh_structure.no_of_modes):
         for vu in range(phh_structure.no_of_modes):
          self.ELPH_sum[iband][k[4]][nu] += np.conjugate(phh_structure.DYN[q_point_no-1][0][mu][nu])*ELPH2[mu][vu]*phh_structure.DYN[q_point_no-1][0][vu][nu] /multik[k[4]]     
  '''




  '''
  for iband in range(self.fermi_nbnd_el):
#   sumweight=0
   for numk, k in enumerate(self.KPOINTS):
      q=self.KPQ[numk]
      ELPH2=np.zeros(shape=(phh_structure.no_of_modes,phh_structure.no_of_modes), dtype=complex) #stores elph for given iband, k [summed over jband
      for jband in range(self.fermi_nbnd_el):
       weight=el_structure.w0gauss(-ell_structure.ENE_fs[jband][q[2]]) *self.WEIGHTS_OF_K[q[1]]
       for iipert in range(phh_structure.no_of_modes):
        for jjpert in range(phh_structure.no_of_modes):
         ELPH2[iipert][jjpert]+=np.conjugate(ELPH[iband][jband][numk][iipert])*ELPH[iband][jband][numk][jjpert]*weight
      ELPH2=ph_structure.symmetrize(phh_structure.nat,np.array(phh_structure.PATT[q_point_no-1]),
ELPH2, structure.at,structure.e, structure.SYMM_crystal, self.SYMM_q)
      for nu in range(phh_structure.no_of_modes):
       for mu in range(phh_structure.no_of_modes):
        for vu in range(phh_structure.no_of_modes):
          self.ELPH_sum[iband][k[4]][nu] += np.conjugate(phh_structure.DYN[q_point_no-1][0][mu][nu])*ELPH2[mu][vu]*phh_structure.DYN[q_point_no-1][0][vu][nu] /multik[k[4]]     

  self.ELPH_sum=self.ELPH_sum #/sumw
  '''

  '''


  multik=np.zeros((len(structure.NONEQ)))
  for k in self.KPOINTS: multik[k[4]]+=1
  for iband in range(self.fermi_nbnd_el):
#   sumweight=0
   for numk, k in enumerate(self.KPOINTS):
   # weight2=el_structure.w0gauss(-ell_structure.ENE_fs[iband][k[4]])  #*self.WEIGHTS_OF_K[k[3]]
    for jband in range(self.fermi_nbnd_el):
     weight=el_structure.w0gauss(-ell_structure.ENE_fs[jband][self.KPQ[numk][2]]) *self.WEIGHTS_OF_K[self.KPQ[numk][1]]
     for nu in range(phh_structure.no_of_modes):
         self.ELPH_sum[iband][k[4]][nu] += np.conjugate(ELPH[iband][jband][numk][nu])*ELPH[iband][jband][numk][nu]*weight /multik[k[4]] #*weight2

  '''



 def elph_single_q_in_whole_kgrid(self,q,structure,phh_structure,ell_structure,l_or_gep):
  print(str(q)+': print elph to file')
  self.lambda_or_elph=l_or_gep
  self.COLORS=np.zeros(shape=(self.fermi_nbnd_el, self.nkp),dtype=complex)
  print(max([i[4] for i in self.KPOINTS]),[len(bnd) for bnd in ell_structure.ENE_fs])
  wagi0=[[el_structure.w0gauss(-ell_structure.ENE_fs[jband][k[4]]) for numk,k in enumerate(self.KPOINTS)] for jband in range(len(self.COLORS)) ]
  wagi=[[wagi0[jband][numk]*self.WEIGHTS_OF_K[numk] for numk in range(len(self.KPOINTS))] for jband in range(len(self.COLORS)) ]
  dos=sum([sum(w) for w in wagi])/sum(self.WEIGHTS_OF_K)
  if self.lambda_or_elph=='elph':
   for numk,k in enumerate(structure.NONEQ):
    for jband in range(self.fermi_nbnd_el): 
     self.COLORS[jband][numk]= np.sum(self.ELPH_sum[jband][numk])
  else:
   for numk,k in enumerate(structure.NONEQ):
    for jband in range(self.fermi_nbnd_el): 
     self.COLORS[jband][numk]+= 2.*(np.sum([(elph)/(phh_structure.FREQ[q-1][phh_structure.ORDER_OF_IRR[q-1][num]])**2*RY_TO_THZ*RY_TO_THZ if phh_structure.FREQ[q-1][phh_structure.ORDER_OF_IRR[q-1][num]]>0 else 0 for num,elph in enumerate(self.ELPH_sum[jband][numk])])) #*wagi0[jband][numk] #/dos  #abs=sqrt(real^2+im^2)
  lambda_sum=sum([sum([ self.COLORS[jband][numk]*self.WEIGHTS_OF_K[numk] for numk in range(len(self.COLORS[jband]))]) for jband in range(len(self.COLORS)) ])
  print (q,': lambda_sum=',lambda_sum,'dos=',dos)


  self.COLORS=[[ bnd[k[3]] for k in structure.allk] for bnd in self.COLORS]

  print (q,":",len(self.COLORS),[len(i) for i in self.COLORS])
  print (q,":",len(ell_structure.ENE_fs),[len(i) for i in ell_structure.ENE_fs])
  print (q,":",len(structure.allk))

def cond(ind,nk):
 if ind[0]!=-1 and ind[0]!=nk[0] and ind[1]!=-1 and ind[1]!=nk[1] and ind[2]!=-1 and ind[2]!=nk[2]: return True
 else: return False

def interpolatee(nk,allk_crystal,ENE):
 h=open('fs/vfermi.frmsf','r')
 tmp=[[float(m) for m in i.split()] for i in h.readlines()]
 h.close()
 nki=[ int(m) for m in tmp[0]]
 nkit=int(nk[0]*nk[1]*nk[2])
 nb=int(tmp[2][0])
 ENEi=np.reshape(tmp[6:6+nkt*nb],(nb,nkt))
 weights=[  [] for i in range(nki[2])  for j in range(nki[1])  for k in range(nki[0])]
 allk_crystal2=[[[ None for i in range(nki[2])] for j in range(nki[1])] for k in range(nki[0])]
 for ni,i in enumerate(allk_crystal):
  allk_crystal2[int(i[0]*nki[0])][int(i[1]*nki[1])][int(i[2]*nki[2])]=ni
 print(ENEi[0])
 pm=[-1,1]
 for i in range(nki[0]):
  for j in range(nki[1]):
   for k in range(nki[2]):
    if allk_crystal2[i][j][k]!=None: weights[i*nki[2]*nki[1]+j*nki[2]+k].append([allk_crystal2[i][j][k],1])
    else:
     for i2 in pm:
      for j2 in pm:
       for k2 in pm
        ind=[i+i2,j+j2,k+k2]
        if cond(ind,nk) and allk_crystal2[ind[0]][ind[1]][ind[2]]!=None: weights[i*nki[2]*nki[1]+j*nki[2]+k].append([allk_crystal2[ind[0]][ind[1]][ind[2]],1])
        for m in weights[i*nki[2]*nki[1]+j*nki[2]+k]: m[1]=m[1]/len(weights[i*nki[2]*nki[1]+j*nki[2]+k])
 return weights,ENEi
    
