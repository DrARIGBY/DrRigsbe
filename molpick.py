import numpy as np
import math 
import random as rn
import scipy as spy
import sys
import scipy.linalg.blas
from Winmol import veiwmol
from CM import CentM
from Analysis_Tools import Center_Of_Mas
__author__='Adam Rigby'
__copyright__ = 'Copyright 2017, Adam M. Rigby'
__credits__ = 'Adam M. Rigby'
__license__ = 'GPL v3.0'
__version__ = '0.1'
__maintainer__ = 'Adam M. Rigby'
__email__ = 'adam.rigby.ar@googlemail.com'
__status__ = 'Development'
def Count_Word(w):
	return sum(ord(c)-64 for c in w)
	
def pickpoints(array1,resi,name,size,box,options,type):
	#~ print 'there are 3 axis needed to calculate the order parameters'
	#~ print 'the long axis, the medium axis and the short axis'
	#~ print 'the long and medium axis will be used to generate the short axis'
	#~ print 'the short and long axis will then be used to generate a new medium axis'
	#~ print 'this is done to make sure that all the axis are othogonal'
	mols=int(len(resi))
	text=''
	#~ print 'mols',mols
	
	a=[[0 for x in range(0,2)] for y in range(0,mols)]
	b=[[0 for x in range(0,2)] for y in range(0,mols)]
	val1=0
	array=np.zeros((mols))
	selected=[0]*4
	Center,Uvecs,MasterVecs,All_Arrays,Total_mas,COM=[],[],[],[],[],[]
	Total_mas,Master_COM=[],[]
	for j in range(0,mols):
		val1=len(array1[j])/size[j]
		array2=spliter(array1[j],size[j],box)
		array=np.asarray(array2)
		#~ print 'splitshape',array.shape
		if options[1]==True:
			selected=veiwmol(array[0],size[j],resi[j],name[j])
			a[j][0],a[j][1],b[j][0],b[j][1]=selected[0],selected[1],selected[2],selected[3]
			Vec=Extract_Vectors(array,a[j],b[j],size[j])
			#~ np.stack((Uvecs[j],Vec))
			#~ print 'vec',Vec
			#~ print 'vec[0]',Vec[0] #this gives the first molecules a, b and c axis vectors
			#~ print 'vec[0][0]',Vec[0][0] #this gives the first molecues a vector
			#~ print 'vec[0][0][0]',Vec[0][0][0] #this gives the first molecules a vectors x value
			Uvecs.append(Vec)
			# #~ print 'Uvecs j',np.asarray(Uvecs[j]).shape # this is the shape of the array of molecule type j
			# #~ print 'Uvecs all',np.asarray(Uvecs[j][0]).shape # this is the shape of j's molcule components (so the number of j molecules and then the 3x3 array for each)
			# #~ print 'Uvecs all',np.asarray(Uvecs[j][0][0]).shape # 3x3 array for each molecule. NOTE there is a Uvecs[j][k][a][v] for a axis vector of type [j] molecule [k],
			# #~ axis eliment [a], vector [v]. And there is a Uvecs[j][k][a][v][x] for the x value of vector[v], axis[a], molecule[k],
			# #~ molecule type [j]. If you np.asarray(Uvecs).shape the whole array, this will create an unusual output.
			
		elif options[1]==False:
			Uvec,Vec=[0],[0]
		#~ print 'name',len(name)
		#~ print 'name[0]',len(name[j] )#this gives the first molecules a, b and c axis vectors
		#~ print 'name[0][0]',len(name[j][0]) #this gives the first molecues a vector
		#~ print 'name[0][0][0]',len(name[j][0][0]) #this give
		#~ print 'array',array
		#~ print 'array[0]',len(array[j] )#this gives the first molecules a, b and c axis vectors
		#~ print 'array[0][0]',len(array[j][0]) #this gives the first molecues a vector
		#~ print 'array[0][0][0]',len(array[j][0][0]) #this give
		Center,mas=Center_Of_Mas(array,name[j],resi[j],options)  
		Total_mas.append(mas)										
		All_Arrays.append(array)	
		COM.append(Center)
	for i in range(0,mols):
		Master_COM.extend(COM[i])
		if options[1]==True:
			MasterVecs.extend(Uvecs[i])
		elif options[1]==False:
			MasterVecs=[0]
	COM=np.asarray(COM)
	All_Arrays=np.asarray(All_Arrays)
	MasterVecs=np.asarray(MasterVecs)
	Uvecs=np.asarray(Uvecs)
	#~ Center=np.asarray(Center)
	return All_Arrays, Uvecs, MasterVecs, COM,Total_mas,Master_COM
	
def spliter(array,size,box):
	array=np.asarray(array)
	val1=len(array)/size #number of molecles 
	#~ print 'val n size',val1,size
	newarray=np.zeros((val1,size,3))
	l=-1
	j=0
	y=len(array)
	#~ print 'val n size',val1,size,y
	for i in range(0,y):
		if i%size==0:
			l+=1
			j=0
			#~ print l,j,i
		
		#~ print l,j,array[i]
		for k in range(0,3):
			newarray[l][j][k]=array[i][k]
		j+=1
	
	newarray=pbc(newarray,size,box)
	#~ newarray=np.asarray(newarray)
	return newarray
	
def pbc(array,size,box):   # ##need to test this properly/THIS DOES NOT CURRENTLY WORK
	title=''
	title+='PBCtest'+str(size)+'.xyz'
	output1=open(title, 'w')
	atomtot=0
	val1=len(array) #molecule number 
	val2=len(array[0]) # #~ number of lenght 
	val3=len(array[0][0]) # #~ xyz components 
	print 'pbc vals',val1,val2
	check,testing=np.zeros((3)),np.zeros((3))
	atomtot=val1*val2
	output1.write('%s \n' % atomtot)
	output1.write('pbc test structure \n')
	#~ check,testing=[[],[],[]],[[],[],[]]
	#~ array2=[]
	array2=np.zeros((val1,val2,3))
	for i in range(0,val1):
		check=array[i][0]
		#~ print 'check',check
		#~ array2.append([])
		for j in range(0,val2):
			for k in range(0,3):
				if (check[k]-array[i][j][k])>box[k]/2.0:
					array2[j][i][k]=array[i][j][k]+box[k]
				elif (check[k]-array[i][j][k])<(-box[k]/2.0):
					array2[i][j][k]=array[i][j][k]-box[k]
				else:
					array2[i][j][k]=array[i][j][k]
			output1.write('C %8.5f  %8.5f  %8.5f \n' % (array2[i][j][0]*10.0, array2[i][j][1]*10.0, array2[i][j][2]*10.0))
			#~ array2[i].append(testing)
			#~ array2[i][j][k]=testing[k]
			#~ print 'testing',testing
			#~ print 'array2',array2[i][j],i,j
		#~ print 'array2',array2[i]
			#~ wait()
		#~ print 'box',box
		#~ print 'array2',array2[i]
	output1.close()
	#~ wait() #YOU NEED TO TEST THIS FIRST
	return array2


def Extract_Vectors(array1,a,b,size):
	aa,bb=np.zeros((2,3)),np.zeros((2,3))
	veca,vecb,vecc=np.zeros((3)),np.zeros((3)),np.zeros((3))
	dya=np.zeros((3,3))
	#~ print suma
	Iden=np.zeros((3,3))
	#~ nomols=int(len(array1))/size
	Iden=np.identity(3)
	n=0
	
	val1=len(array1)
	vectors=np.zeros((val1,3,3))
	for i in range(0,val1): #runs for val1 number of molecules 
		n=0
		for j in range(0,size): #each molecules is of size "size"
			#~ if n==size:
				#~ n=0
			
			if n==a[0]:
				aa[0][0],aa[0][1],aa[0][2]=array1[i][j][0],array1[i][j][1],array1[i][j][2]
				#~ print 'a1[i][j]',array1[i][j],n,j, a[0]
			if n==a[1]:
				aa[1][0],aa[1][1],aa[1][2]=array1[i][j][0],array1[i][j][1],array1[i][j][2]
				#~ print 'a2[i][j]',array1[i][j],n,j, a[1]
			if n==b[0]:
				bb[0][0],bb[0][1],bb[0][2]=array1[i][j][0],array1[i][j][1],array1[i][j][2]
				#~ print 'b1[i][j]',array1[i][j],n,j, b[0]
			if n==b[1]:
				bb[1][0],bb[1][1],bb[1][2]=array1[i][j][0],array1[i][j][1],array1[i][j][2]
				#~ print 'b2[i][j]',array1[i][j],n,j,size, b[1]
			n+=1
		veca=unitvec(aa[0],aa[1])
		#~ print 'veca',veca
		vecb=unitvec(bb[0],bb[1])
		#~ print 'vecb',vecb
		vectors[i][0]=veca
		vecc=np.cross(veca,vecb) #delc is perpendicular to dela and delb
		vecc=Norm_Vec(vecc)
		vectors[i][2]=vecc
		vectors[i][1]=np.cross(vecc,veca) #now delb is perpendicular to dela and delc
		
	return vectors
		
def unitvec(a,b):
	delc=np.zeros((3))
	for i in range(0,3):
		#~ print '(a[i]-b[i])',(a[i]-b[i]),i
		delc[i]=(a[i]-b[i])*(a[i]-b[i])
	hyp=(delc[0]+delc[1]+delc[2])**0.5
	#~ print 'hyp'hyp
	for i in range(0,3):
		delc[i]=(a[i]-b[i])/hyp
	#~ print 'vector out',delc
	return delc 
	
def Norm_Vec(a):
	hyp=(a[0]*a[0]+a[1]*a[1]+a[2]*a[2])**0.5
	return a/hyp
	
	
def Cross_Product(A,B):
	C=np.zeros((3))
	C[0]=A[1]*B[2]-A[2]*B[1]
	C[1]=A[2]*B[0]-A[0]*B[2]
	C[2]=A[0]*B[1]-A[1]*B[0]
	return C
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		