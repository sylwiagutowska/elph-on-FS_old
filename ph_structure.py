import xml.etree.ElementTree as ET
import numpy as np
PRECIS=6



class ph_structure():
 def __init__(self,structure):
  self.FREQ=[] #[i][j] i: q-point; j-nbnd
  self.NONDEG=[]
  self.NONDEG_no=[]
  self.DEG=[] #DEG[q][nbnd] = no of band, with whom the nbnd is degenerated
  self.ORDER_OF_IRR=[]
  self.nat=0
  self.Q=[]
  self.no_of_modes=0
  self.DYN=[]
  self.DYN2=[]
  self.multiplicity_of_qs=[]
#  self.SYMMQ=[]
  self.PATT=[]
  self.prefix=structure.prefix
  self.elph_dir=structure.tmp_dir+'/_ph0/'+self.prefix+'.phsave/'  

 def read_dyn_of_q(self,tmp):
  '''
  for q in range(1,len(self.Q)+1):
     tree = ET.parse(self.elph_dir+'/dynmat.'+str(q)+'.0.xml')
     root = tree.getroot()
     M=( [float(m) for m  in root.find('PARTIAL_MATRIX/PARTIAL_DYN').text.replace(',',' ').split()])
     M2=[]
     for n in range(int(len(M)/2)):
      M2.append(M[2*n]+1j*M[2*n+1])
     DYNMAT=np.reshape(M2,(3,3))
     EIG=np.linalg.eig(DYNMAT)
     self.DYN.append(DYNMAT)
     self.ORDER_OF_IRR.append(np.argsort([abs(i) for i in EIG[0]]))
#     print  (EIG) #([np.sqrt(abs(i)) for i in EIG[0]])
#  exit()
  '''
  all_dyn=[]
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
#     break #if 'Diagonalizing' in i: break  
     all_dyn.append(dyn)
    if 'Diagonalizing' in i: break  
  EIG=np.linalg.eig((all_dyn[0]))
  print('eigenvalue',EIG[0],'multiplicity',len(all_dyn))
  self.DYN.append(all_dyn)


 def read_dyns_nonsymmetrized(self):
  for q in range(1,len(self.Q)+1):
    DYNMATS=np.zeros((3*self.nat,3*self.nat),dtype=complex)
    for i in range(self.no_of_modes+1):
     try: tree = ET.parse(self.elph_dir+'/dynmat.'+str(q)+'.'+str(i)+'.xml')
     except: continue
     root = tree.getroot()
     M=( [float(m) for m  in root.find('PARTIAL_MATRIX/PARTIAL_DYN').text.replace(',',' ').split()])
     M2=[]
     for n in range(int(len(M)/2)):
      M2.append(M[2*n]+1j*M[2*n+1])
     DYNMAT=np.reshape(M2,(3*self.nat,3*self.nat))
     DYNMATS+=DYNMAT
    self.DYN2.append(DYNMATS)  

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
   self.NONDEG_no.append([])
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
   for numf2,f2 in enumerate(self.NONDEG[-1]):
    self.NONDEG_no[-1].append([])
    for nmode,mode in enumerate(self.DEG[-1]):
     if mode==numf2:
      self.NONDEG_no[-1][-1].append(nmode)     
  print (self.NONDEG_no) 
  print (self.FREQ)

  for q in range(1,len(self.Q)+1):
     '''
    DYNMAT=np.zeros((3,3),dtype=complex)
    for i in range(self.no_of_modes):
     try: tree = ET.parse(self.elph_dir+'/dynmat.'+str(q)+'.'+str(i)+'.xml')
     except: continue
     root = tree.getroot()
     M=( [float(m) for m  in root.find('PARTIAL_MATRIX/PARTIAL_DYN').text.replace(',',' ').split()])
     M2=[]
     for n in range(int(len(M)/2)):
      M2.append(M[2*n]+1j*M[2*n+1])
     DYNMAT=np.reshape(M2,(3,3))
     EIG=np.linalg.eig((DYNMAT))
     '''
     EIG=np.linalg.eig((self.DYN[q-1][0]))
#     self.ORDER_OF_IRR.append(np.argsort([abs(i) for i in EIG[0]]))
     self.ORDER_OF_IRR.append([i for i in range(len(EIG[0]))])
#  exit()

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
  self.read_dyns_nonsymmetrized()


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


##################################

def dyn_pattern_to_cart(nat, u, dyn):
  phi=np.zeros((nat,nat,3,3),dtype=complex)
  for i in range(3*nat):
   na=int(np.floor(i/3))
   icart=int(i-3*na)
   for j in range(3*nat):
    nb=int(np.floor(j/3))
    jcart=int(j-3*nb)
    work=0j
    for mu in range(3*nat):
     for nu in range(3*nat):
      work+=u[mu][i] * dyn[nu][mu] * np.conjugate(u[nu][j])
    phi[na][nb][jcart][icart]=work
  return phi

def trn_to_cart(wrk,bg):
  phi=np.zeros((3,3),dtype=complex)
  for i in range(3):
   for j in range(3):
    for k in range(3):
     for l in range(3):
      phi[j][i]+= wrk[l][k]*bg[k][i]*bg[l][j]
  return phi

def trn_to_cryst(phi,at):
  wrk=np.zeros((3,3),dtype=complex)
  for i in range(3):
   for j in range(3):
    for k in range(3):
     for l in range(3):
      wrk[j][i]+= phi[l][k]*at[i][k]*at[j][l]
  return wrk

def compact_dyn(phi,nat):
  dyn=np.zeros((3*nat,3*nat),dtype=complex)
  for na in range(nat):
   for icart in range(3):
    imode=3*na+icart
    for nb in range(nat):
     for jcart in range(3):
      jmode=3*nb+jcart
      dyn[jmode][imode]=phi[na][nb][jcart][icart]
  return dyn

def scompact_dyn(dyn,nat):
  phi=np.zeros((nat,nat,3,3),dtype=complex)
  for na in range(nat):
   for icart in range(3):
    imode=3*na+icart
    for nb in range(nat):
     for jcart in range(3):
      jmode=3*nb+jcart
      phi[na][nb][icart][jcart]=dyn[imode][jmode]
  return phi   

def symmetrize_if_minus_q(phi,s,nt,irotmq):
 for na in range(nat):
  for nb in range(nat): 
   for ipol in range(3):
    for jpol in range(3):
     work=np.zeros((3,3),dtype=complex)
     sna=na #irt[na][irotmq]
     snb=nb #irt[nb][irotmq]
     arg=0
     for kpol in range(3):
      arg+= 0 #q[kpol] * (rtau[na][irotmq][kpol] - rtau[nb][irotmq][kpol])
     arg=arg*np.pi
     fase=np.cos(arg)+1j*np.sin(arg)
     for kpol in range(3):
      for lpol in range(3):
       work[jpol][ipol] += s[irotmq][kpol][ipol]*s[irotmq][lpol][jpol]* phi[sna][snb][lpol][kpol]*fase
     phip[na][nb][jpol][ipol]=(phi[na][nb][jpol][ipol]+ np.conjugate(work[jpol][ipol]))*0.5
 return phip

def symmetrize_small_qgroup(phi,s,nat):
 sinv=[np.linalg.inv(i) for i in s]
 iflb=np.zeros((nat,nat),dtype=complex)
 faseq=np.zeros((len(s)),dtype=complex)
 for na in range(nat):
  for nb in range(nat): 
   if iflb[na][nb]==0: 
    work=np.zeros((3,3),dtype=complex)
    for isymq in range(len(s)):
     irot=isymq
     sna=na #irt[na][irot]
     snb=nb #irt[nb][irot]
     arg=0 #sum([ q[ipol]* (rtau[na][irot][ipol] - rtau[nb][irot][ipol]) for ipol in range(3)])
     arg=arg*np.pi
     faseq[isymq]=1 #np.cos(arg)+1j*np.sin(arg)
     for ipol in range(3):
      for jpol in range(3):
       for kpol in range(3):
        for lpol in range(3):
         work[jpol][ipol] += s[irot][kpol][ipol]*s[irot][lpol][jpol] *phi[sna][snb][lpol][kpol] * faseq[isymq]
    for isymq in range(len(s)):
      irot=isymq
      sna=na #irt[na][irot]
      snb=nb #irt[nb][irot]
      for ipol in range(3):
       for jpol in range(3):
        phi[sna][snb][jpol][ipol]=0.j
        for kpol in range(3):
         for lpol in range(3):
          phi[sna][snb][jpol][ipol] += sinv[irot][kpol][ipol]*sinv[irot][lpol][jpol] * work[lpol][kpol]*np.conjugate(faseq[isymq])
    iflb[sna][snb]=1
 phi=phi/len(s)
 return phi  


def symmetrize(nat,pattern,dyn,at,bg,s,s_of_q ):
  phi=dyn_pattern_to_cart(nat, pattern, dyn)
  for na in range(nat):
   for nb in range(nat):
    phi[na][nb]=trn_to_cryst(phi[na][nb],at)
  #impose hermecity
  for na in range(nat):
   for nb in range(nat):
    for ipol in range(3):
     for jpol in range(3):
      phi[na][nb][ipol][jpol]=0.5 * (phi[na][nb][ipol][jpol] + np.conjugate(phi[nb][na][jpol][ipol]))
      phi[nb][na][jpol][ipol]=np.conjugate(phi[na][nb][ipol][jpol])
  phi=symmetrize_small_qgroup(phi,s_of_q,nat)
  for na in range(nat):
   for nb in range(nat):
    phi[na][nb]=trn_to_cart(phi[na][nb],bg)
  phi=compact_dyn(phi,nat)
  return phi.transpose()


