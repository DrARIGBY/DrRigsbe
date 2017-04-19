import numpy as np
from numpy import linalg as LA
from scipy import linalg as SLA
import math
__author__='Adam Rigby'
__copyright__ = 'Copyright 2017, Adam M. Rigby'
__credits__ = 'Adam M. Rigby'
__license__ = 'GPL v3.0'
__version__ = '0.1'
__maintainer__ = 'Adam M. Rigby'
__email__ = 'adam.rigby.ar@googlemail.com'
__status__ = 'Development'
def Center_Of_Mas(array,atname,resi,option):
	cent1=[]
	if option[0] == True:
	#~ while option:
		text1='CENTMas'+str(resi).replace('[','').replace(']','').replace("'" ,'')+'.xyz'
		text2='CENTNo'+str(resi).replace('[','').replace(']','').replace("'" ,'')+'.xyz'
		text3='CENTR'+str(resi).replace('[','').replace(']','').replace("'" ,'')+'.xyz'
		output1=open(text1, 'w')
		output2=open(text2, 'w')
		output3=open(text3, 'w')
		y=np.genfromtxt("atomlist.txt",dtype=None,usecols=(0,1,2))
		cent1,cent2,name,mass,label,rads=[],[],[],[],[],[]
		val1=len(array) # number of molecles 
		val2=len(array[0]) # molecule length
		print'val2',val2
		print 'vals cmas loop',val1,val2
		#~ print 'sghape of atname',atname
		output1.write('%s \n' % val1)
		output1.write('Center of Mass calculated with mass structure \n')
		output2.write('%s \n' % val1)
		output2.write('Center of Mass calculated without mass structure \n')
		output3.write('%s \n' % val1)
		output3.write('Center of Mass with diameter + VDW radius in angstrom \n')
		for i in range(len(y)):
			mass.append(y[i][0])
			name.append(y[i][1])
			label.append(y[i][2])
		cent1,cent2=np.zeros((val1,3)),np.zeros((val1,3))
		rads,hyp=np.zeros((val1)),np.zeros((3))
		totmas=0.0
		weight=0
		#~ print'len atname',len(atname)
		#~ print'len atname[0]',len(atname[0])
		#~ print'len atname[0][0]',len(atname[0][0])
		#~ wait()
		for i in range(0,val2):
			for k in range(0,len(y)): #number of atoms masses to compare
				#~ print 'atname[i]',atname[i],i,k,len(y)
				
				check2=str(atname[i]).replace('[','').replace(']','').replace("'",'')
				check1=str(label[k]).replace('[','').replace(']','').replace("'",'')
				#~ print'checks',check1,check2
				if check2 == check1:
					weight=k
					totmas+=mass[weight]
		for i in range(0,val1): #no. of molecules
			
			pos1,pos2=np.zeros((3)),np.zeros((3))
			#~ print'i',i
			#~ cent.append([])
			#~ print 'at len',atname.shape
			for j in range(0,val2): #no of atoms per molecule
				
				for k in range(0,3):
					pos1[k]+= array[i][j][k]*mass[weight]
					pos2[k]+= array[i][j][k]
				for l in range(j,val2):
					for k in range(0,3):
						hyp[k]=(array[i][j][k]-array[i][l][k])
					radius=((hyp[0]**2.0+hyp[1]**2.0+hyp[2]**2.0)**0.5)/2.0
					if radius>rads[i]:
						rads[i]=radius
		#~ print 'radius of i',rads[i],i
			for k in range(0,3):
				cent1[i][k]=(pos1[k]/totmas)
				cent2[i][k]=(pos2[k]/val2)
			output1.write('C %8.5f  %8.5f  %8.5f \n' % (cent1[i][0]*10.0, cent1[i][1]*10.0, cent1[i][2]*10.0))
			output2.write('C %8.5f  %8.5f  %8.5f \n' % (cent2[i][0]*10.0, cent2[i][1]*10.0, cent2[i][2]*10.0))
			output3.write('C %8.5f  %8.5f  %8.5f  %8.5f   \n' % (cent1[i][0]*10.0, cent1[i][1]*10.0, cent1[i][2]*10.0,rads[i]*10.0+1.7))
		#~ print 'out of the cm loop'
	#~ wait()
	else:
		val1=len(array) # number of molecles 
		val2=len(array[0]) # molecule length
		#~ print 'vals cmas loop',val1,val2
		cent1=np.zeros((val1,3))
		rads,hyp=np.zeros((val1)),np.zeros((3))
		for i in range(0,val1): #no. of molecules
			pos2=np.zeros((3))
			for j in range(0,val2): #no of atoms per molecule
				for k in range(0,3):
					pos2[k]+= array[i][j][k]
				for l in range(j,val2):
					for k in range(0,3):
						hyp[k]=(array[i][j][k]-array[i][l][k])
					radius=((hyp[0]**2.0+hyp[1]**2.0+hyp[2]**2.0)**0.5)/2.0
					if radius>rads[i]:
						rads[i]=radius
			#~ print 'radius of i',rads[i],i
			for k in range(0,3):
				cent1[i][k]=(pos2[k]/val2)
			#~ print 'cent1',cent1[i]
		#~ print 'cent1 in COM',np.asarray(cent1).shape
		totmas=val1
	return  np.asarray(cent1),totmas

def Order_Parameters(vectors): #extract the unitvectors and do the dyadic product	
	EigVal,selectvec=np.zeros((3,3)),np.zeros((3,3))
	Eig,sum_Q,Qtens,dya=np.zeros((3,3,3)),np.zeros((3,3,3)),np.zeros((3,3,3)),np.zeros((3,3,3))
	Iden,select=np.identity(3),np.zeros((3))
	Qx,Qy,Qz=np.zeros((3,3)),np.zeros((3,3)),np.zeros((3,3))
	n=0
	Points_Selected=[0]*3
	vectors=np.asarray(vectors)
	val1=len(vectors)
	for i in range(0,val1):
		#~ print ''
		for j in range(0,3):
			dya[j]=np.outer(vectors[i][j],vectors[i][j])
		for j in range(0,3):
			sum_Q[j]=sum_Q[j]+(3.0*dya[j]-Iden)
	for i in range(0,3):
		#~ print ''
		#~ print 'NEW Q CALCULATION',i
		Qtens[i]=sum_Q[i]/(2.0*val1)  # convert these to matrix
		EigVal[i],Eig[i]=LA.eigh(Qtens[i],'U') #Generates Eigen values (EigVal) and Eigen vectors (Eig)
		select[i]=np.amax(EigVal[i]) #Eigen values are checked and the max value is selected for each set of EigVals
		for j in range(0,3):
			if select[i]==EigVal[i][j] :
				#~ selectvec[i]=Eig[i][j]
				selectvec[i][0]=-Eig[i][0][j]
				selectvec[i][1]=-Eig[i][1][j]
				selectvec[i][2]=-Eig[i][2][j]
	MaxEig=np.amax(select) #selects the largest of the 3 largest eignevlas
	MinEig=np.amin(select)
	for i in range(0,3):
		if  select[i] != MaxEig and select[i] !=MinEig :
			MidEig=select[i]
			Points_Selected[1]=i
			Qy=Qtens[i]
		if select[i]==MaxEig :
			Points_Selected[0]=i
			Qz=Qtens[i]
		if select[i]==MinEig :
			Points_Selected[2]=i
			Qx=Qtens[i]

	ZQzZ=np.dot(selectvec[Points_Selected[0]],(np.dot(Qz,selectvec[Points_Selected[0]].T)))
	YQyY=np.dot(selectvec[Points_Selected[1]],(np.dot(Qy,selectvec[Points_Selected[1]].T)))
	XQxX=np.dot(selectvec[Points_Selected[2]],(np.dot(Qx,selectvec[Points_Selected[2]].T)))
	YQxY=np.dot(selectvec[Points_Selected[1]],(np.dot(Qx,selectvec[Points_Selected[1]].T)))
	XQyX=np.dot(selectvec[Points_Selected[2]],(np.dot(Qy,selectvec[Points_Selected[2]].T)))
	Q_00=ZQzZ
	Q_22=(XQxX+YQyY-YQxY-XQyX)/3

	Zarray=np.zeros((val1))
	print 'Q_00',ZQzZ
	print 'Q_22',YQyY
	for i in range(0,val1):
		Zarray[i]=math.degrees(Angle_Between(selectvec[Points_Selected[0]],vectors[i][Points_Selected[0]]))
	return Zarray

def Radial_Distribution_Function(COM,box,resi,mas):
	names=len(resi)
	label=''
	for i in range(0,names):
		label+=resi[i]
	text1='RDF'+str(label).replace('[','').replace(']','').replace("'" ,'')+'.xyz'
	#~ print 'COM',COM
	output1=open(text1, 'w')
	print 'Generating the RDF file called',text1
	box=np.asarray(box)
	maxrad=(box[0]**2.0+box[1]**2.0+box[2]**2.0)**0.5
	bin_size=(maxrad/500.0)
	radial_bin=[0]*500
	dr=np.zeros((3))
	seg=np.zeros((3))
	val1=len(COM)
	volume=box[0]*box[1]*box[2]
	density=val1*mas/(box[0]*box[1]*box[2])
	mastot=val1*mas
	print 'density',density,'bin_size',bin_size
	for i in range(0,val1-1):
		#~ print 'COM[i]',COM[i]
		for j in range(i+1,val1):
			#~ print 'COM[j]',COM[j]
			#~ if j>i:
			for k in range(0,3):
				seg[k]=((COM[i][k]-COM[j][k])**2.0)**0.5
				dr[k]=(box[k]*(math.fabs(seg[k]/box[k]-round(seg[k]/box[k],0))))
				#~ if seg[k] > box[k]/2:
					#~ dr[k] = seg[k]-box[k]
				#~ else:
					#~ dr[k]=seg[k]
				#~ dr[k]=dr[k]+box[k]*(math.fabs(seg[k]/box[k]))
			#~ print 'dr', dr 
			#~ print 'COM I',COM[i]
			#~ print 'COM J',COM[j]
			#~ print 'radial_value 1',round((((dr[0]**2.0+dr[1]**2.0+dr[2]**2.0)**0.5)/bin_size),0)
			radial_value=int(round((((dr[0]**2.0+dr[1]**2.0+dr[2]**2.0)**0.5)/bin_size),0))
			#~ print 'radial_value used',radial_value
			radial_bin[radial_value]+=2
	for i in range(0,500):
		if i*bin_size<box[0]/2.0:
			if radial_bin[i] !=0:
				radial_bin[i]=(radial_bin[i])/(density*4.0*(np.pi)*(((i+1)*bin_size)**2.0)*bin_size)
			else:
				radial_bin[i]=0
			output1.write('%8.5f  %8.5f \n' % (i*bin_size, radial_bin[i]))

	
def Print_COMS_Buddies(array,pairs,buddies,resi,type,types_groups,paired_groups):
	label=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn']
	text1='AllCOM'+str(type)
	numtot=0
	for i in range(0,len(resi)):
		text1+=str(resi[i]).replace('[','').replace(']','').replace("'",'')
		numtot+=int(len(array[i]))
	#~ print 'numtot',numtot
	text1+='.xyz'
	output1=open(text1, 'w')
	output1.write(str(numtot)+'\n')
	output1.write(str(type)+' Structure with '+ str(text1)+'\n')
	for i in range(0,int(len(array))):
		for j in range(0,int(len(array[i]))):
			loop=0
			for p in range(0,pairs):
				for n in range(0,buddies+1): 
					if (types_groups[i][p][n]==i) and (paired_groups[i][p][n]==j) and loop==0:
						loop=p
						#~ print'p',p
					elif loop==0:
						loop=0
			output1.write('%3s %8.5f  %8.5f  %8.5f \n' % (label[loop], array[i][j][0]*10.0, array[i][j][1]*10.0, array[i][j][2]*10.0))
	

def History_Radial_Distribution(COM,box,resi,pairs,budies,rsmall_all,rsmall_groups,rsmall_same,paired_all,paired_groups,paired_same,types_all,types_groups,types_same):
	pass
	
	
def Parallel_RDF(COM,box,resi,mas):
	names=len(resi)
	label=''
	for i in range(0,names):
		label+=resi[i]
	box=np.asarray(box)
	radial_bin=[0]*500
	segment=20
	segments=20.0
	x_split=box[0]/segments
	y_split=box[1]/segments
	z_split=box[2]/segments
	maxrad=(x_split**2.0+box[1]**2.0+box[2]**2.0)**0.5
	volume=x_split*box[1]*box[2]
	bin_size=(maxrad/500.0)
	val1=len(COM)
	#~ print 'val1',val1
	#~ wait()
	COM=np.asarray(COM)
	#~ print 'com',COM
	#~ print 'com',COM.shape
	dr=np.zeros((3))
	seg=np.zeros((3))
	for i in range(0,segment):
		delx1=x_split*i
		delx2=x_split*(i+1)
		
		found,density=0,0
		for j in range(0,val1):
			if COM[j][0]>delx1 and COM[j][0]<delx2:
				for k in range(0,val1):
					#~ print 'CMj',COM[j][0],'CMk',COM[k][0]
					if k!=j and COM[k][0]>delx1 and COM[k][0]<delx2:
						#~ print 'dels',delx1,delx2
						#~ print 'CMj',COM[j][0],'CMk',COM[k][0]
						for l in range(0,3):
							seg[l]=((COM[k][l]-COM[j][l])**2.0)**0.5
							dr[l]=(box[l]*(math.fabs(seg[l]/box[l]-round(seg[l]/box[l],0))))
						radial_value=int(round((((dr[0]**2.0+dr[1]**2.0+dr[2]**2.0)**0.5)/bin_size),0))
						radial_bin[radial_value]+=2
						#~ print 'radial_bin',radial_bin[radial_value], radial_value
						found+= 1
		if found!=0:
			#~ print 'density',(found*mas/volume),i
			text1='PXRDF'+str(label).replace('[','').replace(']','').replace("'" ,'')+ str(i+1) +'.xyz'
			print 'Generating parallel RDF file called',text1
			output1=open(text1, 'w')
			for j in range(0,500):
				if j*bin_size<box[0]/2.0:
					radial_bin[j]=radial_bin[j]/((val1/volume)*4.0*(np.pi)*(((j+1)*bin_size)**2.0)*bin_size)
					output1.write('%8.5f  %8.5f \n' % (j*bin_size, radial_bin[j]*segments))
			output1.close()
	
	radial_bin=[0]*500	
	for i in range(0,segment):
		delx1=y_split*i
		delx2=y_split*(i+1)
		
		found,density=0,0
		for j in range(0,val1):
			if COM[j][1]>delx1 and COM[j][1]<delx2:
				for k in range(0,val1):
					if k!=j and COM[k][1]>delx1 and COM[k][1]<delx2:
						#~ print 'dels',delx1,delx2
						#~ print 'CMj',COM[j][0],'CMk',COM[k][0]
						for l in range(0,3):
							seg[l]=((COM[k][l]-COM[j][l])**2.0)**0.5
							dr[l]=(box[l]*(math.fabs(seg[l]/box[l]-round(seg[l]/box[l],0))))
						radial_value=int(round((((dr[0]**2.0+dr[1]**2.0+dr[2]**2.0)**0.5)/bin_size),0))
						radial_bin[radial_value]+=2
						#~ print 'radial_bin',radial_bin[radial_value], radial_value
						found+= 1
		if found!=0:
			#~ print 'density',(found*mas/volume),i
			text1='PYRDF'+str(label).replace('[','').replace(']','').replace("'" ,'')+ str(i+1) +'.xyz'
			print 'Generating parallel RDF file called',text1
			output1=open(text1, 'w')
			for j in range(0,500):
				if j*bin_size<box[0]/2.0:
					radial_bin[j]=radial_bin[j]/((val1/volume)*4.0*(np.pi)*(((j+1)*bin_size)**2.0)*bin_size)
					output1.write('%8.5f  %8.5f \n' % (j*bin_size, radial_bin[j]*segments))
			output1.close()
			
	radial_bin=[0]*500
	for i in range(0,segment):
		delx1=z_split*i
		delx2=z_split*(i+1)
		found,density=0,0
		for j in range(0,val1):
			if COM[j][2]>delx1 and COM[j][2]<delx2:
				for k in range(0,val1):
					if k!=j and COM[k][2]>delx1 and COM[k][2]<delx2:
						#~ print 'dels',delx1,delx2
						#~ print 'CMj',COM[j][0],'CMk',COM[k][0]
						for l in range(0,3):
							seg[l]=((COM[k][l]-COM[j][l])**2.0)**0.5
							dr[l]=(box[l]*(math.fabs(seg[l]/box[l]-round(seg[l]/box[l],0))))
						radial_value=int(round((((dr[0]**2.0+dr[1]**2.0+dr[2]**2.0)**0.5)/bin_size),0))
						radial_bin[radial_value]+=2
						#~ print 'radial_bin',radial_bin[radial_value], radial_value
						found+= 1
		if found!=0:
			#~ print 'density',(found*mas/volume),i
			text1='PZRDF'+str(label).replace('[','').replace(']','').replace("'" ,'')+ str(i+1) +'.xyz'
			print 'Generating parallel RDF file called',text1
			output1=open(text1, 'w')
			for j in range(0,500):
				if j*bin_size<box[0]/2.0:
					radial_bin[j]=radial_bin[j]/((val1/volume)*4.0*(np.pi)*(((j+1)*bin_size)**2.0)*bin_size)
					output1.write('%8.5f  %8.5f \n' % (j*bin_size, radial_bin[j]*segments))
			output1.close()
	
def Dot_Product_Qii(a,Q,b):
	Q_ii=b[0]*(a[0]*Q[0,0]+a[1]*Q[1][0]+a[2]*Q[2][1])+b[1]*(a[0]*Q[0][1]+a[1]*Q[1][1]+a[2]*Q[2][1])+b[2]*(a[0]*Q[0][2]+a[1]*Q[1][2]+a[2]*Q[2][2])
	return Q_ii

def Dyadic_Product_Qii(A,B):
	Q=np.zeros((3,3))
	for i in range(0,3):
		for j in range(0,3):
			Q[i][j]=A[i]*B[j]
	return Q

def Angle_Between(v1,v2):
	v1u=Unit_Vec(v1)
	v2u=Unit_Vec(v2)
	ang=np.arccos(np.dot(v1u,v2u))
	#~ print 'ang in angle between',ang,np.dot(v1u,v2u)
	#~ print 'v1 v2', v1u,v2u
	#~ print 'from angle_between', np.arccos(np.clip(np.dot(v1u,v2u),-1.0,1.0))
	return ang

def Unit_Vec(vector):
	return vector / np.linalg.norm(vector)

	#~ paired=np.zeros((types,pairs,budies,2))
	#~ print'1len(paired)',len(paired)
	#~ print'1len(paired[0])',len(paired[0])
	#~ print'1len(paired[0][0])',len(paired[0][0])
	#~ print'1len(paired[0][0][0])',len(paired[0][0][0])
	#~ print'1len(paired[0][0][0][0])',len(paired[0][0][0][0])
	
	#~ paired=[[[[0]*2]*budies]*pairs]*types
	#~ print'2len(paired)',len(paired)
	#~ print'2len(paired[0])',len(paired[0])
	#~ print'2len(paired[0][0])',len(paired[0][0])
	#~ print'2len(paired[0][0][0])',len(paired[0][0][0])
	#~ print'2len(paired[0][0][0][0])',len(paired[0][0][0][0])
	
	#~ paired=[[[[0 for x in range(2)] for y in range(budies)] for z in range(pairs)] for zz in range(types)]
	#~ print'3len(paired)',len(paired)
	#~ print'3len(paired[0])',len(paired[0])
	#~ print'3len(paired[0][0])',len(paired[0][0])
	#~ print'3len(paired[0][0][0])',len(paired[0][0][0])
	#~ print'3len(paired[0][0][0][0])',len(paired[0][0][0][0])
	#~ pare_label=[[0]*pairs]*types 

	