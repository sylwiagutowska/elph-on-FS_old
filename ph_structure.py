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
  self.Q_crystal=[]
  self.no_of_modes=0
  self.DYN=[]
  self.DYN2=[]
  self.multiplicity_of_qs=[]
  self.qstar=[]
#  self.SYMMQ=[]
  self.PATT=[]
  self.prefix=structure.prefix
  self.prefixdyn=structure.prefixdyn
  self.elph_dir=structure.tmp_dir+'/_ph0/'+self.prefix+'.phsave/'  
  self.e=structure.e

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
         dyn[mu][nu]=tmp_dyn[na*self.nat+nb][ii][jj] #/sqrt(mass_a*mass_b)/amu_ry
#     break #if 'Diagonalizing' in i: break  
     all_dyn.append(dyn)
    if 'Diagonalizing' in i: break  
 # EIG=np.linalg.eig((all_dyn[0]))
#  print('eigenvalue',EIG[0],'multiplicity',len(all_dyn))
  self.DYN.append(all_dyn)


 def read_dyns_nonsymmetrized(self):
  for q in range(1,len(self.Q)+1):
    DYNMATS=np.zeros((3*self.nat,3*self.nat),dtype=complex)
    for i in range(3*self.nat+1):
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
     self.Q.append(np.round([float(m) for m in i.split()[3:6] ],PRECIS))
     qcryst=np.round(np.dot(np.linalg.inv(np.transpose(self.e)),self.Q[-1]),PRECIS)
     for m in range(3):
      while qcryst[m]<0: qcryst[m]+=1
      while qcryst[m]>=1: qcryst[m]-=1
     self.Q_crystal.append(np.round(qcryst,PRECIS-1))
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
 # print (self.NONDEG_no) 
 # print (self.FREQ)

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
  #   self.ORDER_OF_IRR.append(np.argsort([abs(i) for i in EIG[0]]))
     self.ORDER_OF_IRR.append([i for i in range(len(EIG[0]))])
#  exit()

 def read_multiplicity_of_q(self,dynfile):
  multiplicity=0
  qstars=[]
  for ni,i in enumerate(dynfile):
   if 'Dynamical  Matrix in cartesian axes' in i:
    multiplicity=multiplicity+1
    qstars.append([float(m) for m in dynfile[ni+2].split()[3:6]])
  self.qstar.append(qstars)
  self.multiplicity_of_qs.append( multiplicity)

 def read_ph_structure(self):
  ###read frequencies and make DEG[q][nbnd] = no of band, with whom the nbnd is degenerated
  print(' read phonon structure')
  h=open(self.prefixdyn+'.dyn1')
  for i in range(3): tmp=h.readline()
  h.close()
  self.nat=int(tmp.split()[1])
  for file in range(1,100):
   try: h=open(self.prefixdyn+'.dyn'+str(file))
   except: break
   tmp=h.readlines()
   h.close()
   self.read_freq_of_q(tmp)
   self.read_dyn_of_q(tmp)
   self.read_q(tmp)
   self.read_multiplicity_of_q(tmp)
  self.check_degeneration_of_modes()
  self.read_dyns_nonsymmetrized()



 def read_patterns(self):
  for q in range(len(self.Q)):
   self.PATT.append([])
   tree = ET.parse(self.elph_dir+'/patterns.'+str(q+1)+'.xml')
   root = tree.getroot()
   nirr=int(root.find('IRREPS_INFO/NUMBER_IRR_REP').text)
   for i in range(nirr):
     rep=root.find('IRREPS_INFO/REPRESENTION.'+str(i+1))
     npert=int(rep.find('NUMBER_OF_PERTURBATIONS').text)
     for j in range(npert):
      pat=[ (lambda m: complex(float(m[0]),float(m[1])))(n.replace(',',' ').split()) for n in  rep.find('PERTURBATION.'+str(j+1)+'/DISPLACEMENT_PATTERN').text.split('\n')[1:-1]]
      self.PATT[-1].append(pat)


##################################

def dyn_pattern_to_cart(nat, u, dyn):
#  ! This routine rotates the dynamical matrix dynwork written
#  ! in cartesian basis to the basis of the patterns u and adds it to
#  ! the dynamical matrix dyn that is supposed to be in the basis of the
#  ! patterns.
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
      work+= u[mu][i] * dyn[nu][mu] * np.conjugate(u[nu][j])
    phi[nb][na][jcart][icart]=work
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
      dyn[jmode][imode]=phi[nb][na][jcart][icart]
  return dyn

def scompact_dyn(dyn,nat):
  phi=np.zeros((nat,nat,3,3),dtype=complex)
  for na in range(nat):
   for icart in range(3):
    imode=3*na+icart
    for nb in range(nat):
     for jcart in range(3):
      jmode=3*nb+jcart
      phi[nb][na][jcart][icart]=dyn[jmode][imode]
  return phi   

def symmetrize_if_minus_q(phi,s,nat,irt,irotmq,rtau,q):
 phip=np.zeros(phi.shape,dtype=complex)
 for na in range(nat):
  for nb in range(nat): 
   for ipol in range(3):
    for jpol in range(3):
     work=np.zeros((3,3),dtype=complex)
     sna=irt[na][irotmq]
     snb=irt[nb][irotmq]
     dr=rtau[na][irotmq] - rtau[nb][irotmq]
     arg=sum([ q[kpol]* dr[kpol] for kpol in range(3)])
     arg=arg*2*np.pi
     fase=np.cos(arg)+1j*np.sin(arg)
     for kpol in range(3):
      for lpol in range(3):
       work[jpol][ipol] += s[irotmq][kpol][ipol]*s[irotmq][lpol][jpol]* phi[snb][sna][lpol][kpol]*fase
     phip[nb][na][jpol][ipol]=(phi[nb][na][jpol][ipol]+ np.conjugate(work[jpol][ipol]))*0.5
 return phip

def symmetrize_small_qgroup(phi,s,nat,irt,rtau,q):
 sinv=[np.linalg.inv(i) for i in s]
 iflb=np.zeros((nat,nat),dtype=complex)
 faseq=np.zeros((len(s)),dtype=complex)
 for na in range(nat):
  for nb in range(nat): 
   if iflb[nb][na]==0: 
    work=np.zeros((3,3),dtype=complex)
    for isymq in range(len(s)):
     irot=isymq
     sna=irt[na][irot]
     snb=irt[nb][irot]
     dr=rtau[na][irot] - rtau[nb][irot]
     arg=sum([ q[ipol]* dr[ipol] for ipol in range(3)])
     arg=arg*2*np.pi
     faseq[isymq]=np.cos(arg)+1j*np.sin(arg)
     for ipol in range(3):
      for jpol in range(3):
       for kpol in range(3):
        for lpol in range(3):
         work[jpol][ipol] += s[irot][kpol][ipol]*s[irot][lpol][jpol] *(phi[snb][sna][kpol][lpol]) * faseq[isymq]
    for isymq in range(len(s)):
      irot=isymq
      sna=irt[na][irot]
      snb=irt[nb][irot]
      for ipol in range(3):
       for jpol in range(3):
        phi[snb][sna][jpol][ipol]=0.j
        for kpol in range(3):
         for lpol in range(3):
          phi[snb][sna][jpol][ipol] += sinv[irot][kpol][ipol]*sinv[irot][lpol][jpol] * (work[lpol][kpol])*np.conjugate(faseq[isymq])
      iflb[snb][sna]=1
 phi=phi/len(s)
 return phi  

def if_minus_q(q,s):
 mq=-q
# for i in range(3):
#  while mq[i]<0: mq[i]+=1
#  while mq[i]>=1: mq[i]-=1
 for ns,sym in enumerate(s):
  mq2=np.round(np.dot(sym,mq),PRECIS)
  for i in range(3):
   while mq2[i]<0: mq2[i]+=1
   while mq2[i]>=1: mq2[i]-=1
  mq2=np.round(mq2,PRECIS)
  q2=np.round(q,PRECIS)
  if mq2[0]==q2[0] and mq2[1]==q2[1] and mq2[2]==q2[2]: return ns
 else: return -1

def impose_hermicity(phi,nat):
  #impose hermecity
  for na in range(nat):
   for nb in range(nat):
    for ipol in range(3):
     for jpol in range(3):
      phi[nb][na][jpol][ipol]=0.5 * (phi[nb][na][jpol][ipol] + np.conjugate(phi[na][nb][ipol][jpol]))
      phi[na][nb][ipol][jpol]=np.conjugate(phi[nb][na][jpol][ipol])
  return phi 

def symmetrize(nat,pattern,dyn,at0,bg,s,s_of_q,irt,rtau,q,q_cart ):
  #rtau= input: the R associated at each  ,  irt (48, nat)=input: the rotated of each atom, isym=s_of_q=input: the small group of q
  at=np.linalg.inv(np.transpose(bg)) #/2/np.pi
  phi=dyn_pattern_to_cart(nat, pattern, dyn)
  for na in range(nat):
   for nb in range(nat):
    phi[nb][na]=trn_to_cryst(phi[nb][na],at)
  phi=impose_hermicity(phi,nat)
  irotmq=if_minus_q(q,s_of_q)
  if irotmq!=-1: phi=symmetrize_if_minus_q(phi,s_of_q,nat,irt,irotmq,rtau,q ) 
  phi=symmetrize_small_qgroup(phi,s_of_q,nat,irt,rtau,q )
  if irotmq!=-1: print(q,'MQ present')
  for na in range(nat):
   for nb in range(nat):
    phi[nb][na]=trn_to_cart(phi[nb][na],bg)
 #   print(na,nb,'\n',phi[na][nb])
  phi=compact_dyn(phi,nat)
  return phi #.transpose()


