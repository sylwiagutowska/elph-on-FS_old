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

 def make_kpoints_single_q(self,file,basic_structure):
  tree = ET.parse(self.elph_dir+'elph.'+str(file)+'.1.xml')
  root = tree.getroot()
  self.nkp=int(root.find('PARTIAL_EL_PHON/NUMBER_OF_K').text)
  self.KPOINTS=[ [ round(float(m),PRECIS)\
         for m in root.find(
          'PARTIAL_EL_PHON/K_POINT.'+str(k)+'/COORDINATES_XK'
          ).text.split() ]\
           for k in range(1,self.nkp+1)]
  for i in range(self.nkp): self.KPOINTS[i].append(i)
  structure_new=structure.structure()
  structure_new.NONEQ=self.KPOINTS
  structure_new.SYMM=basic_structure.SYMM
  structure_new.e=basic_structure.e
  structure_new.no_of_kpoints=basic_structure.no_of_kpoints
  structure_new.make_kgrid()
  self.KPOINTS_all=structure_new.allk


 def read_elph_single_q(self,file,ph_structure,el_structure): #read info from first file
  tree = ET.parse(self.elph_dir+'elph.'+str(file)+'.1.xml')
  root = tree.getroot()
  self.nbnd_el=int(root.find('PARTIAL_EL_PHON/NUMBER_OF_BANDS').text)
  ELPH=[[[ [] for k in range(self.nkp)] for j in range(self.nbnd_el)] for i in range(self.nbnd_el)]
  self.ELPH_sum=[[ [] for k in range(self.nkp)] for j in range(self.nbnd_el) ]

  #read elph and dynmat from all files
  for mode in range(1,len(ph_structure.NONDEG[file-1])+1):
   #elph
   tree = ET.parse(self.elph_dir+'elph.'+str(file)+'.'+str(mode)+'.xml')
   root = tree.getroot()
   for country in root.iter('PARTIAL_EL_PHON'):
    for k in range(1,self.nkp+1):
     for town in country.iter('K_POINT.'+str(k)):
      partial_elph=town.find('PARTIAL_ELPH')
      if k==1: npert=int(partial_elph.get('size'))/self.nbnd_el/self.nbnd_el
      elph_k=[ [float(m.replace(',',' ').split()[0]), 
                float(m.replace(',',' ').split()[1])] 
                for m in town.find('PARTIAL_ELPH').text.split('\n') if len(m.split())>0  ]
      for jband in range(self.nbnd_el):
       for iband in range(self.nbnd_el):
        for iipert in range(npert):
         ELPH[jband][iband][k-1].append\
          (elph_k[jband*self.nbnd_el*npert+iband*npert+iipert])

      for jband in range(self.nbnd_el):
       for iipert in range(npert):
        self.ELPH_sum[jband][k-1].append([0 for skl in range(2)])
        for iband in range(self.nbnd_el):
         for skl in range(2):
          self.ELPH_sum[jband][k-1][-1][skl]+=          (elph_k[jband*self.nbnd_el*npert+iband*npert+iipert][skl])

 def elph_single_q_in_whole_kgrid(self,q,structure,ph_structure,el_structure):
 #choose only bands which cross EF and sum over j in pairs <i,j>
  print('From all '+str(self.nbnd_el)+' bands detected in elph calc. only bands ',el_structure.bands_num,' cross EF and will be written in frmsf')
  COLORS=[ [] for jband in el_structure.bands_num]
  for numk,k in enumerate(self.KPOINTS):
   for numjband,jband in enumerate(el_structure.bands_num): 
    COLORS[numjband].append( sum([elph[0]**2+elph[1]**2 for elph in self.ELPH_sum[jband][numk]])   )
  
  print len(COLORS),[len(i) for i in COLORS]
  print len(el_structure.ENE),[len(i) for i in el_structure.ENE]
  print len(structure.allk)
  h=open('elph'+str(q)+'.frmsf','w')
  h.write(str(structure.no_of_kpoints[0])+' '+str(structure.no_of_kpoints[1])+' ' +str(structure.no_of_kpoints[2])+'\n')
  h.write('1\n'+str(len(el_structure.bands_num))+'\n')
  for i in structure.e:
   for j in i:
    h.write(str(j)+' ')
   h.write('\n')
  for bnd in el_structure.ENE:
   for k in structure.allk:
    h.write(str(bnd[k[3]])+'\n')
  for bnd in COLORS:
   for k in self.KPOINTS_all:
    h.write(str(bnd[k[3]])+'\n')
  h.close()


