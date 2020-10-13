
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


el_structure=el_structure.el_structure()
el_structure.read_el_structure()


ph_structure=ph_structure.ph_structure()
ph_structure.read_ph_structure()


elph_structure=elph_structure.elph_structure()
elph_structure.make_kpoints_single_q(2,structure)
elph_structure.read_elph_single_q(2,ph_structure,el_structure)
elph_structure.elph_single_q_in_whole_kgrid(2,structure,ph_structure,el_structure)
###
