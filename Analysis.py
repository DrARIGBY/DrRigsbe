from __future__ import print_function
from Input_Files import Gromacs_Gro as GI
from Input_Files import DL_Poly_Input
from Input_Files import XYZ_Big as XYZI
from Input_Files import QMGA_Input as QMGA
from Input_Files import PDB_Input as PDBI
from Input_Files import XYZ_Converter
from Analysis_Tools import Paired_Radial_Distribution as PRD
from Analysis_Tools import Buddy_Radial_Distribution as BRD
from Analysis_Tools import History_Radial_Distribution as HRD
from System import OutPut_Radial_Distribution_All_Buddies as ORDAB
import argparse
import sys, os
import numpy as np
import math 
import random as rn
import scipy as spy
from molpick import pickpoints,Count_Word
from Analysis_Tools import Order_Parameters as OP
from gromacsgro import Gromacs_Gro
from Analysis_Tools import Radial_Distribution_Function,Parallel_RDF

__author__='Adam Rigby'
__copyright__ = 'Copyright 2017, Adam M. Rigby'
__credits__ = 'Adam M. Rigby'
__license__ = 'GPL v3.0'
__version__ = '0.1'
__maintainer__ = 'Adam M. Rigby'
__email__ = 'adam.rigby.ar@googlemail.com'
__status__ = 'Development'
#~ from __future__ import *
if __name__=='__main__':
	"""
	Reads in any input type (pdb,xyz,gro) and generates different setup structures.
	"""
#~ print sys.argv[0]
#~ print sys.argv
	parser = argparse.ArgumentParser(description='Reads a lattice PDB and two molecular PDBs. Replaces lattice points with molecules in either a '
	'pure, blended, or mixed system. Assumes that the interfacial plane for a bileyer is normal to the '
	'Z-axis')
	parser.add_argument('-Input',nargs=1,help='the input file needed to do all the analysis',required=True)
	parser.add_argument('-Analysis',nargs='*',help='name the analsysis methods you wish to use',
		choices=['RDF','COM','com','Order Parameters','Parallel RDF','all','prdf','Center of Mas','Radial Distribution Function','rdf','his','buds'],required=True)
	parser.add_argument('-mols',nargs='*',help='number of molecule')
	parser.add_argument('-num',nargs='*',help='number of total molecules (specificly for QMGA files)')
	parser.add_argument('-Hist-tpr',nargs=1,help='the gromacs file.tpr to go with the file.tpr.trr file')
	parser.add_argument('-Hist-list',nargs=2,help='the gromacs file(x).gro (with out the number-x) and the number of files being used (assumed starting from 0-x)')
	parser.add_argument('-Budds',nargs=1,help='Number of closest molecules to the center of mas of the closest pairs of molecules')
	parser.add_argument('-Pairs',nargs=1,help='Number of closest molecule pairs')
	
	
	args=parser.parse_args()
	
	sets=len(vars(args)['Analysis'])
	velocity,name,atname,pos,resi,box=[],[],[],[],[],[]
	Each_Vector,Master_Vector, Master_Arrays, COM,Master_COM=[],[],[],[],[]
	Zarray=[]
	#~ options=[0]*
	#~ if sets > 0:
	#~ for i in range(0,sets):
	listed_options=['com','op','rdf','prdf','Order Parameters','Parallel RDF','Center of Mas','Radial Distribution Function','his','buds','all']
	name_options=['Center of Mas','Order Parameters','Radial Distribution Function','Parallel Radial Distribution Functions','History RDFs','Paired RDF though History','Buddies']
	options=[False]*len(listed_options)
	print(len(listed_options))
	filename=vars(args)['Input'][0]
	text=(vars(args)['Input'][0]).split('.')
	texts_lists=len(text)
	file_type=''
	for i in range(0,texts_lists):
		text[i]=text[i].replace(' ','').replace('[','').replace(']','').replace("'",'')
		if text[i]=='gro':
			file_type=text[i]
			#~ n+=1
		elif text[i]=='xyz':
			file_type=text[i]
			#~ n+=1
		elif text[i]=='pdb':
			file_type=text[i]
			#~ n+=1
		elif text[i]=='qmga':
			file_type=text[i]
					
	if file_type== '':
		file_type='dlpoly'
	
	#~ if options_list>2:
	print('You have selected these options')
	for i in range(0,sets):
		#~ print 'i',i
		amount1=Count_Word(vars(args)['Analysis'][i].lower())
		#~ print 'amount1',amount1
		for j in range(0,len(listed_options)):
			amount2=Count_Word(listed_options[j].lower())
			#~ print 'amount2',amount2
			if vars(args)['Analysis'][i].lower() in listed_options[-1]:
				options=[True]*len(listed_options)
				for k in range(0,(len(listed_options)-1)):
					print ('Generating',name_options[k])
				break
			elif  vars(args)['Analysis'][i].lower() in listed_options[j] and amount1==amount2 :
				options[j]=True
				if j==6 or j == 0:
					l=0
				elif j==1 or j == 4:
					l=1
				elif j==2 or j == 7:
					l=2
				elif j==3 or j == 5:
					l=3
				elif j==8:
					l=4
				elif j == 9:
					l=5
				print ('Generating',name_options[l])
				
	#~ option
	#~ wait()
	#~ print 'filename' ,filename
	try:
		mol_types=vars(args)['mols'][0]
	except(ValueError, RuntimeError,TypeError,NameError,IOError):
		mol_types=1
	
	if file_type=='gro':
		resi,name, position, velocity, box, size= GI(filename)
		#~ resi,name, position, size=Gromacs_Input(filename)
		print('gro files')
	elif file_type=='pdb':
		resi,name, position, velocity, box, size=PDBI(filename)
		print('pdb files')
	elif file_type=='xyz':
		resi,name,position,size,box=XYZI(filename,mol_types)
		print('xyz files')
	elif file_type=='dlpoly':
		resi,name,position,size,box=DL_Poly_Input(filename,mol_types)
	elif file_type=='qmga':
		try:
			mol_types=vars(args)['num'][0]
		except(ValueError, RuntimeError,TypeError,NameError,IOError):
			mol_types=729
		resi,name,position,size,box=QMGA(filename,mol_types)
	#~ resi,name, position, velocity, box, size= Gromacs_Gro(filename)
	#~ print 'strings1',strings1
	#~ print 'strings2',strings2
	#~ print('options',name)
	#~ print'box',box
	#~ print('position',position.shape)
	molecule_types=len(resi)
	Master_Arrays, Each_Vector , Master_Vector, COM, mas, Master_COM=pickpoints(position,resi,name,size,box,options,file_type)
	#~ print('COM',COM)
	print ('')
	print ('IN MAIN ANALYSIS LOOP')
	XYZ_Converter(position,name,resi)
	#~ print ('Master_Arrays',np.asarray(Master_Arrays).shape)
	#~ print ('master_vector',np.asarray(Master_Vector).shape	)
	#~ print ('each_vector',np.asarray(Each_Vector).shape)
	#~ print ('each_vector[0]',np.asarray(Each_Vector[0]).shape)
	#~ print ('COM',np.asarray(COM).shape)
	#~ print ('COM[0]',np.asarray(COM[1]).shape)
	#~ print ('Master_COM',np.asarray(Master_COM).shape	)
	#~ print ('Master_COM[0]',np.asarray(Master_COM[0]).shape)

	for i in range(0,molecule_types):
		#~ print ''
		#~ print 'ENTERING THE OP LOOP'
		if options[1]==True or options[4]==True  :
			array=OP(Each_Vector[i])
		if options[2]==True or options[7]==True :
			Radial_Distribution_Function(COM[i],box,resi[i],mas[i])
		if options[3]==True or options[5]==True :
			Parallel_RDF(COM[i],box,resi[i],mas[i])
	
	if options[9]==True or options[10]==True :
		try:
			budies=vars(args)['Budds'][0]
		except(ValueError, RuntimeError,TypeError,NameError,IOError):
			budies=3
		try:
			pairs=vars(args)['Pairs'][0]
		except(ValueError, RuntimeError,TypeError,NameError,IOError):
			pairs=5
		
		rsmall_all,rsmall_groups,rsmall_same,paired_all,paired_groups,paired_same,types_all,types_groups,types_same,all_pos,groups_pos,same_pos=PRD(COM,box,resi,pairs,budies,Master_Arrays,size)
		#~ print('types_all main',types_all)
		#~ print('paired_all main',paired_all)
		#~ print('rsmall_all main',rsmall_all)
	if options[9]==True:
		rsmall_all,rsmall_groups,rsmall_same,paired_all,paired_groups,paired_same,types_all,types_groups,types_same=BRD(COM,box,resi,pairs,budies,rsmall_all,rsmall_groups,rsmall_same,paired_all,paired_groups,paired_same,types_all,types_groups,types_same,Master_Arrays,all_pos,groups_pos,same_pos,size)
		#~ print('types_all budies loop',types_all)
		#~ print('paired_all budies loop',paired_all)
		#~ print('rsmall_all budies loop',rsmall_all)
		ORDAB(Master_Arrays,COM,resi,rsmall_groups,paired_groups,types_groups,size,box,budies,pairs)
	if options[10]==True:
		rsmall_all,rsmall_groups,rsmall_same,paired_all,paired_groups,paired_same,types_all,types_groups,types_same=HRD(COM,box,resi,pairs,budies,rsmall_all,rsmall_groups,rsmall_same,paired_all,paired_groups,rsmall_same,types_all,types_groups,rsmall_same)
		
#~ def XYZ_Converter(pos,name,resi):
	#~ for i in range(0,len(resi)):
		#~ text1='Converted_'+str(resi[i]).replace('[','').replace(']','').replace("'" ,'')+ str(i+1) +'.xyz'
		#~ print ('Converting origianl file to xyz under the name: \n',text1)
		#~ output1=open(text1, 'w')
		#~ for j in range(0,len(pos[i])):
			#~ output1.write('%4d %8.5f %8.5f  %8.5f \n' % (name[i][j], pos[i][j][0], pos[i][j][1], pos[i][j][2]))
	
#~ gmx trjconv -f nptdata.tpr.trr -s nptdata.tpr  -t0 0 -e 500 -dt 1 -o test.gro -sep 

		
		
		
		
		
		
		
		