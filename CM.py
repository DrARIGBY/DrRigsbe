import numpy as np
import math
__author__='Adam Rigby'
__copyright__ = 'Copyright 2017, Adam M. Rigby'
__credits__ = 'Adam M. Rigby'
__license__ = 'GPL v3.0'
__version__ = '0.1'
__maintainer__ = 'Adam M. Rigby'
__email__ = 'adam.rigby.ar@googlemail.com'
__status__ = 'Development'
def CentM(position,centa=[]):
	#~ for i in range(len(position)):
		
	tot=[]
	i=0
	#~ for i in range(step):
	y=len(position)
	#~ print y
	tot=np.sum(position, axis=0)
	
	#~ print tot
	centa=np.array([tot/y])
	i=i+1
	#~ print centa
	return centa

def Den(density,molnum):
	volume=molnum/density
	perline=math.pow(molnum, (1.0/3.0))
	linelenght=math.pow(volume, (1.0/3.0))
	return perline, linelenght 


def wtsplit(molnum,totalnum, mol=[], num=[]): #need to create a cheaking loop.
	
	#~ mol.append([])
	#~ num.append([])
	error=0.01
	del1=totalnum*error
	del2=totalnum-del1
	del3=totalnum+del1
	while True:
		mols=[]
		nums=[]
		#~ for i in range(0,molnum):
		nums, mols=split(molnum,totalnum,nums,mols)
			#~ num.extend(num1)
			#~ mols.extend(mol1)
		print 'dels', del1, del2, del3
		print 'initial values'
		print nums
		print mols
		wtper=sum(mols)
		number=int(sum(nums))
		print 'wtper', wtper
		print 'number', number
		if ( wtper != 1):
			print 'Im afraid the weight % dose not sum to 100%'
			del mols[:]
			del nums[:]
			#~ mol=[]
			#~ num=[]
			#~ mol.append([])
			#~ num.append([])
			print 'you will need to re-enter the percentages.'
			#~ print mol
			#~ print 'num 2'
			#~ print num
		else:
			False
			break
	for i in range(0,molnum):
		num.append(nums[i])
		mol.append(mols[i])
	print 'your system contains a mixture with the percentages'
	print mol
	print 'This is equal to this number of moleucles for each type (wtsplit)'
	print num	
		
	return num, mol
			
	
def split(molnum,totalnum,num=[],mol=[]):
	for j in range(0,molnum):
		molper=float(input("what is the weight percentage of molecule %d: \n " % (j+1)))
		#~ num=[]
		#~ mol=[]
		#~ mol.append([])
		#~ num.append([])
		print 'molper',molper
		mol.append(molper)
		print 'mol',mol
		check=totalnum*molper % 1.0
		if check > 0.5 :
			num.append(math.ceil(totalnum*molper))
		else :
			num.append(math.ceil(totalnum*molper)-1)
		print 'in split loop', num
	return num, mol



