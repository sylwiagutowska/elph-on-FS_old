
import xml.etree.ElementTree as ET
import numpy as np
from operator import itemgetter
import structure
import el_structure
import ph_structure
import elph_structure

structure=structure.structure()
structure.read_structure()
structure.make_kgrid()


el_structure=el_structure.el_structure(structure)
el_structure.read_el_structure()


ph_structure=ph_structure.ph_structure(structure)
ph_structure.read_ph_structure()
#ph_structure.check_symm_of_q(structure)
ph_structure.read_patterns()



elph_structure=elph_structure.elph_structure(ph_structure)
print ph_structure.Q
for q in range(1,len(ph_structure.Q)+1):
 print('calculations for '+str(q)+'. of total '+str(len(ph_structure.Q))+' q points')
 elph_structure.make_kpoints_single_q(q,structure)
 elph_structure.read_elph_single_q(q,ph_structure,el_structure)
 elph_structure.elph_single_q_in_whole_kgrid(q,structure,\
                ph_structure,el_structure,'lambda')  #'lambda' or 'elph'
elph_structure.sum_over_q(ph_structure,structure,el_structure)

###

