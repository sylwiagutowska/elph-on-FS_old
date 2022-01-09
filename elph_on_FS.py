
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
el_structure.read_extra_el_structure()
structure_extra= copy.deepcopy(sstructure)
structure_extra.NONEQ=el_structure.NONEQ_extra
structure_extra.WK=el_structure.WK_extra
structure_extra.no_of_kpoints=el_structure.nk_extra
print(structure_extra.NONEQ[:50])
print(sstructure.NONEQ[:50])
structure_extra.make_kgrid()
#structure_extra.find_newkpoints_in_old_list(sstructure.allk)
#structure_extra.find_weights_in_old_list(sstructure.allk)


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
pb_mass=207.2*9.3975038



for qno in range(len(phh_structure.Q)):
 structure_new=copy.deepcopy(sstructure)
 structure_new.check_symm(phh_structure.Q[qno],sstructure.NONEQ)
#print(len(structure_new.SYMM))
# print(phh_structure.DYN2[qno],np.sum(phh_structure.DYN2[qno],axis=0))
 dyn=ph_structure.symmetrize(phh_structure.nat,np.array(phh_structure.PATT[qno]),
phh_structure.DYN2[qno],
sstructure.at,sstructure.e,
 sstructure.SYMM_crystal, structure_new.SYMM_crystal)
 dyn=dyn/pb_mass/AMU_RY
# print(round_complex(dyn))
 print('diagonalized')
 phh_structure.DYN[qno]=[dyn]
 a=np.linalg.eig(dyn)
 phh_structure.FREQ[qno]=( abs(a[0].real)**0.5-(abs(a[0].imag))**0.5 )*RY_TO_THZ
 print( qno+1,phh_structure.DYN[qno], phh_structure.FREQ[qno])
 #print(round_complex(a[0]*RY_TO_THZ*RY_TO_THZ))
 #print(round_complex(a[1]))
exit()


elphh_structure=elph_structure.elph_structure(phh_structure,'lambda')
print (phh_structure.Q)
print(len(phh_structure.Q),len(phh_structure.DYN2),len(phh_structure.PATT))

'''
for q in range(1,len(ph_structure.Q)+1):
 print('calculations for '+str(q)+'. of total '+str(len(ph_structure.Q))+' q points')
 elph_structure.make_kpoints_single_q(q,structure,ph_structure.Q[q-1])
 elph_structure.read_elph_single_q(q,ph_structure,el_structure,structure)
 elph_structure.elph_single_q_in_whole_kgrid(q,structure,\
                ph_structure,el_structure,'lambda')  #'lambda' or 'elph'
elph_structure.sum_over_q(ph_structure,structure,el_structure)
'''
#elphh_structure.parallel_job(sstructure,structure_extra,el_structure,phh_structure)
elphh_structure.single_job([2,sstructure,structure_extra,phh_structure,el_structure,'lambda'])
elphh_structure.sum_over_q(phh_structure,sstructure,el_structure)
elphh_structure.extrapolate_values(sstructure,el_structure)

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

