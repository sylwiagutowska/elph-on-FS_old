
import xml.etree.ElementTree as ET
import numpy as np
from operator import itemgetter
import structure
import el_structure
import ph_structure
import elph_structure
import copy
def round_complex(y):
 prec=9
 try:
  for i in range(len(y)):
   for j in range(len(y[i])):
    y[i][j]=round(y[i][j].real,prec)+1j*round(y[i][j].imag,prec)
 except:
  try:
   for i in range(len(y)):
    y[i]=round(y[i].real,prec)+1j*round(y[i].imag,prec)
  except:
    y=round(y.real,prec)+1j*round(y.imag,prec)
 return y

sstructure=structure.structure()
sstructure.read_structure()
sstructure.make_kgrid()


el_structure=el_structure.el_structure(sstructure)
el_structure.read_el_structure()


phh_structure=ph_structure.ph_structure(sstructure)
phh_structure.read_ph_structure()
#ph_structure.check_symm_of_q(structure)

phh_structure.read_patterns()
#print(phh_structure.PATT[1])
basic_structure=sstructure


 #check if composing dynmat works - WORKS! GIVES PROPER omega^2

ELECTRONMASS_SI  = 9.10938215e-31   # Kg
AMU_SI           = 1.660538782e-27  #Kg
AMU_AU           = AMU_SI / ELECTRONMASS_SI
AMU_RY           = AMU_AU / 2. #=911.44
RY_TO_THZ=3289.8449

#pb_mass=207.2 #*(9.3975038) #**2
#pbbi_mass=207.2 #*6.6485187 #**2
amass=basic_structure.at_masses #[pbbi_mass,pbbi_mass]
print(amass)
#print(phh_structure.Q_crystal)
#sstructure.calc_irt()
#print(sstructure.irt)
#exit()
jedynki=[]
h=open('macierze_dyn.dat','w')
for qno in range(len(phh_structure.Q)):
 print(qno,phh_structure.Q_crystal[qno])
 structure_new=copy.deepcopy(sstructure)
 structure_new.check_symm(phh_structure.Q[qno],sstructure.NONEQ)
 structure_new.calc_irt()
#print(len(structure_new.SYMM))
# print(phh_structure.DYN2[qno],np.sum(phh_structure.DYN2[qno],axis=0))
 dyn=ph_structure.symmetrize(phh_structure.nat,np.array(phh_structure.PATT[qno]),
     (phh_structure.DYN2[qno]), sstructure.at,sstructure.e,
     sstructure.SYMM_crystal, structure_new.SYMM_crystal,
     structure_new.irt,structure_new.rtau,phh_structure.Q_crystal[qno],phh_structure.Q[qno])
# print(round_complex(dyn))
 print('symmetrizied')
 for i in range(3):
  for na in range(phh_structure.nat):
   mu=3*(na)+i
   for j in range(3):
    for nb in range(phh_structure.nat):
     nu=3*(nb)+j
     dyn [mu][nu] = dyn [mu][nu] / ( amass[nb]*amass[na])**0.5

 dyn=np.round(dyn,decimals=10)
 dyn=dyn/AMU_RY
 phh_structure.DYN[qno]=[dyn]
 a=np.linalg.eigvals(dyn)
# phh_structure.DYN[qno]=[np.transpose(a[1])]
 prevfreq=phh_structure.FREQ[qno]
 phh_structure.FREQ[qno]=np.sqrt(np.abs(a)) *RY_TO_THZ
 for ni,i in enumerate(a):
  if i<0: phh_structure.FREQ[qno][ni]=-phh_structure.FREQ[qno][ni]
 phh_structure.FREQ[qno]=sorted(phh_structure.FREQ[qno])
 if qno==0: 
  for m in range(3): phh_structure.FREQ[qno][m]=0
 h.write(str(qno)) 
 for i in phh_structure.DYN[qno][0]:
  for j in i:
   h.write(str(j)+' ')
  h.write('\n')
 jedynki.append([round(i/prevfreq[ni],4) for ni,i in enumerate(sorted(phh_structure.FREQ[qno]))])
 for i in jedynki[-1]:  h.write(str(i)+'\n')
 if sum(jedynki[-1])!=len(jedynki[-1]): print(structure_new.SYMM_crystal,phh_structure.Q_crystal[qno],phh_structure.Q[qno])
h.close()
 #print(round_complex(a[0]*RY_TO_THZ*RY_TO_THZ))
 #print(round_complex(a[1]))


for ni,i in enumerate(jedynki):
 if sum(i)!=len(i): print('Q',ni+1,phh_structure.Q_crystal[ni],'wrong! freq/freq=',i)

exit()

elphh_structure=elph_structure.elph_structure(phh_structure,'lambda')



'''
for q in range(1,len(ph_structure.Q)+1):
 print('calculations for '+str(q)+'. of total '+str(len(ph_structure.Q))+' q points')
 elph_structure.make_kpoints_single_q(q,structure,ph_structure.Q[q-1])
 elph_structure.read_elph_single_q(q,ph_structure,el_structure,structure)
 elph_structure.elph_single_q_in_whole_kgrid(q,structure,\
                ph_structure,el_structure,'lambda')  #'lambda' or 'elph'
elph_structure.sum_over_q(ph_structure,structure,el_structure)
'''
elphh_structure.parallel_job(sstructure,el_structure,phh_structure)
#elphh_structure.single_job([26,sstructure,phh_structure,el_structure,'lambda'])
elphh_structure.sum_over_q(phh_structure,sstructure,el_structure)
#elphh_structure.extrapolate_values(sstructure,el_structure)

'''
h=open('lambda.frmsf')
tmp=h.readlines()[6+len(el_structure.ENE_fs)*len(structure.allk):]
h.close()
m=0
elph_structure.SUMMED_COLORS=[[ 0 for i in range(len(structure.allk))] for j in range(len(el_structure.ENE_fs))]

for bnd in range(len(el_structure.ENE_fs)):
 for k in range(len(structure.allk)):
    elph_structure.SUMMED_COLORS[bnd][k]=float(tmp[m])
    m+=1
print(m)
elph_structure.extrapolate_values(structure,el_structure)
'''
###

