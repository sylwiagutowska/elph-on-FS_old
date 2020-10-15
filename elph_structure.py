import xml.etree.ElementTree as ET
import numpy as np
import structure
PRECIS=6
class elph_structure():
 def __init__(self):
  self.ELPH_sum=[]
  self.elph_dir='tmp_dir/_ph0/ir.phsave/'
  self.KPOINTS=[]
  self.nbnd_el=0
  self.nkp=0
  self.KPOINTS_all=[]

 def make_kpoints_single_q(self,q,basic_structure):
  self.pm=[0,-1,1]
  tree = ET.parse(self.elph_dir+'elph.'+str(q)+'.1.xml')
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
  

 def w0gauss(self,x):
  degauss=0.02
  sqrtpm1= 1. / 1.77245385090551602729
  # cold smearing  (Marzari-Vanderbilt-DeVita-Payne)
  arg = min (200., (x - 1.0 /  (2.0)**0.5 ) **2.)
  return sqrtpm1 * np.exp ( - arg) * (2.0 - ( 2.0)**0.5 * x)/0.02


 def dyn_to_cart(self,nat,pat,matrix):
  DYN_PHI=[[[[complex(0,0) for i in range(3)] for j in range(3)] for nb in range(3*nat)] for na in range(3*nat)]
  for i in range(nat):
     na =int( i / 3) 
     icart = i - 3 * (na )
     for j in range(3*nat):
        nb = int(j/3)
        jcart = j - 3 * (nb)
        work = complex(0,0)
        for mu in range(3*nat):
         for nu in range(3*nat):
          work = work +pat[mu][i]*matrix[mu][nu] * pat[nu][j].conjugate()
        DYN_PHI [na][nb][icart][jcart] = work
  return DYN_PHI

 def transform_matrix_to_crystal(self,mat,at):
  return [[ sum([ sum([ mat[k][l]*at[k][i]*at[l][j] for l in range(3)]) for k in range(3)]) for j in range(3)] for i in range(3)]

 def transform_matrix_to_cartesian(self,mat,e):
  return [[ sum([ sum([ mat[k][l]*e[i][k]*e[j][l] for l in range(3)]) for k in range(3)]) for j in range(3)] for i in range(3)]

 def symmetrize(phi,nat,minus_q)
  !    We start by imposing hermiticity
  for na in range(nat):
   for nb in range(nat):
    for ipol in range(3):
     for jpol in range(3):
      phi[na][nb][ipol][jpol]= 0.5*(phi[na][nb][ipol][jpol]+(phi[nb][na][jpol][ipol].conjugate() )
      phi[nb][na][jpol][ipol] =phi[na][nb][ipol][jpol].conjugate()
  !    Then we impose the symmetry q -> -q+G if present
  if (minus_q):
   for na in range(nat):
    for nb in range(nat):
     work=[[complex(0,0) for i in range(3)] for j in range(3)]
     for ipol in range(3):
      for jpol in range(3):
                 sna = irt (irotmq, na)
                 snb = irt (irotmq, nb)
                 arg = 0.d0
                 do kpol = 1, 3
                    arg = arg + (xq (kpol) * (rtau (kpol, irotmq, na) - &
                                              rtau (kpol, irotmq, nb) ) )
                 enddo
                 arg = arg * tpi
                 fase = CMPLX(cos (arg), sin (arg) ,kind=DP)
                 do kpol = 1, 3
                    do lpol = 1, 3
                       work (ipol, jpol) = work (ipol, jpol) + &
                            s (ipol, kpol, irotmq) * s (jpol, lpol, irotmq) &
                            * phi (kpol, lpol, sna, snb) * fase
                    enddo
                 enddo
                 phip (ipol, jpol, na, nb) = (phi (ipol, jpol, na, nb) + &
                      CONJG( work (ipol, jpol) ) ) * 0.5d0
              enddo
           enddo
        enddo
     enddo
     phi = phip
  endif


 def read_elph_single_q(self,q,ph_structure,el_structure,structure): #read info from first file
  tree = ET.parse(self.elph_dir+'elph.'+str(q)+'.1.xml')
  root = tree.getroot()
  self.nbnd_el=int(root.find('PARTIAL_EL_PHON/NUMBER_OF_BANDS').text)
  ELPH=[[[ [] for k in range(self.nkp)] for j in range(self.nbnd_el)] for i in range(self.nbnd_el)] #stores elph[k][ibnd][jbnd][nmode]
  self.ELPH_sum=[[ [] for k in range(self.nkp)] for j in range(self.nbnd_el) ]
#stores elph[k][ibnd][nmode][mmode]

  #read elph  from all files
  for mode in range(1,len(ph_structure.NONDEG[q-1])+1):
   #elph
   tree = ET.parse(self.elph_dir+'elph.'+str(q)+'.'+str(mode)+'.xml')
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

  for jband in range(self.nbnd_el):
    for k in range(1,self.nkp+1):
       self.ELPH_sum[jband][k-1]=[[complex(0,0) for i in range(ph_structure.no_of_modes)] for j in range(ph_structure.no_of_modes)]
       for iipert in range(ph_structure.no_of_modes):
        for jjpert in range(ph_structure.no_of_modes):
         for iband in range(self.nbnd_el):
           self.ELPH_sum[jband][k-1][iipert][jjpert]+=\
                self.w0gauss(el_structure.ef-el_structure.ENE[jband][self.KPOINTS[k-1][4]])*\
                ELPH[jband][iband][k-1][iipert].conjugate()*\
                ELPH[jband][iband][k-1][jjpert]
  for jband in range(self.nbnd_el):
    for k in range(1,self.nkp+1):
      self.ELPH_sum[jband][k-1]=\
               self.dyn_to_cart(ph_structure.nat,ph_structure.PATT[q-1],\
                                self.ELPH_sum[jband][k-1])
      self.ELPH_sum[jband][k-1]=\
		self.transform_matrix_to_crystal(self.ELPH_sum[jband][k-1],structure.at)

 def elph_single_q_in_whole_kgrid(self,q,structure,ph_structure,el_structure):
 #choose only bands which cross EF and sum over j in pairs <i,j>
  print('From all '+str(self.nbnd_el)+' bands detected in elph calc. only bands ',el_structure.bands_num,' cross EF and will be written in frmsf')
  COLORS=[ [] for jband in el_structure.bands_num]
  for numk,k in enumerate(self.KPOINTS):
   for numjband,jband in enumerate(el_structure.bands_num): 
    COLORS[numjband].append( sum([elph.real**2+elph.imag**2 for elph in self.ELPH_sum[jband][numk]])   )
  
  print len(COLORS),[len(i) for i in COLORS]
  print len(el_structure.ENE_fs),[len(i) for i in el_structure.ENE_fs]
  print len(structure.allk)
  h=open('elph'+str(q)+'.frmsf','w')
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


