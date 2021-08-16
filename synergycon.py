# -*- coding: utf-8 -*-
"""
Created on Mon May 31 01:29:05 2021

@author: zhang
"""

import sys
import numpy as np
from pyitlib import discrete_random_variable as drv
from numpy import mat;
import itertools
import pandas as pd


def funccdpconsygy(x,y):
	X0=mat(x)
	X1 = np.array(y)
	X1 = [int(numeric_string) for numeric_string in X1]
	[X0rows, X0cols] = X0.shape
	X02=[]
	for line in range(X0rows):
		stre=''
		line = X0[line,:]
		thestate = stre.join(np.array(line).tolist()[0])
		X02.append(thestate)
	X02 = pd.DataFrame(X02)
	X02 = X02[0].astype('category')
	X0 = list(X02.cat.codes)
	Iout = drv.entropy(X1, base=2)-drv.entropy_conditional(X=X1,Y=X0, base=2)
	return(Iout)

def funpvaleconsygy(amat, bmat, y, thisv, iternum=10):
	mycount=0
	cmat = np.column_stack((amat,bmat))
	thresholdv = thisv
	arange = amat.shape[1]
	brange = bmat.shape[1]
	if iternum > 0:
		for xx in range(0, iternum):
			amat = shflmat(amat,arange)
			bmat = shflmat(bmat,brange)
			cmat = np.column_stack((amat,bmat))
			sumres = funccdpconsygy(cmat,y)
			res1 = funccdpconsygy(amat,y)
			res2 = funccdpconsygy(bmat,y)
			thisv = sumres - res1 -res2
			if thisv > thresholdv:
				mycount = mycount+1	   
		pvalue = (mycount+1) / (iternum + 1)
		return(pvalue)
	else:
		return('')

def digitcon(xx,digitnumber):
	xx = np.array(xx)
	xx = [float(numeric_string) for numeric_string in xx]
	xx = np.digitize(xx, np.linspace(min(xx), max(xx), digitnumber))
	xx = [int(numeric_string) for numeric_string in xx]
	return(xx)


def shflmat(x,time):
	for ii in range(time):
		np.random.shuffle(x[:,ii])
	return(x)

def stroutconsygy(x,y,proteinname):
	x = [proteinname[str(numeric_string)] for numeric_string in x]
	y = [proteinname[str(numeric_string)] for numeric_string in y]
	str1 = '_'
	x = str1.join(x)
	y = str1.join(y)
	out = x + ':' + y
	return(out)

def synergycon(datain,outputfile,dataanno, digitnumber=2,iternum=0, limitinter=5):
	ii = 0
	proteinname={}
	with open(dataanno, "r") as f:
		for line in f.readlines():
			line = line.strip('\n')  
			line = line.split("\t")
			proteinname1 = line[1]
			proteinname[str(ii)] = proteinname1
			ii = ii + 1
	X0=[]
	X1=[]
	with open(datain, "r") as f:
		for line in f.readlines():
			line = line.strip('\n')  
			line = line.split("\t")
			distype = line[-1]
			line = line[:-1]
			X0.append(line)
			X1.append(distype)
	
	X0=mat(X0)

	[X0rows, X0cols] = X0.shape
	for line in range(X0cols):
		thestate = np.array(X0[:,line])
		thestate = thestate.flatten()
		thestate = digitcon(thestate,digitnumber)
		X0[:,line] = np.array([thestate]).T
	X0=np.array(X0)
	y=X1

	snpnum = X0.shape[1]
	if snpnum < 3:
		return "too small matrix"
	naoutput = (0,1)
	
	namelist = list(range(snpnum))
	namelist = [[ct] for ct in namelist]
	treemat = np.random.rand(snpnum,snpnum)
	treemat[:] = np.nan
	allcb = list(itertools.combinations(list(range(len(namelist))), 2))
	output = []
	for i2 in allcb:
		fi = namelist[i2[0]]
		ai = namelist[i2[1]]
		alli = fi+ai
		sumres = funccdpconsygy(X0[:,alli],y)
		res1 = funccdpconsygy(X0[:,fi],y)
		res2 = funccdpconsygy(X0[:,ai],y)
		thisv = sumres - res1 -res2
		treemat[fi,ai] = thisv
	inds = np.unravel_index(np.nanargmax(treemat, axis=None), treemat.shape)
	thev = treemat[inds]
	if thev < 0:
		output = mat(output)
		np.savetxt(outputfile,output,fmt='%s',delimiter='\t',newline='\n')
		return(1)
	thep = funpvaleconsygy(X0[:, namelist[inds[0]] ], X0[:, namelist[inds[1]] ], y, thev, iternum)
	out = stroutconsygy(namelist[inds[0]],namelist[inds[1]],proteinname)
	out = [out, thev, thep]
	output.append(out)
	namelist[inds[0]]= namelist[inds[0]]+namelist.pop(inds[1])
	lenmerged = len(namelist[inds[0]])
	treemat = np.delete(treemat,inds[1],axis=0)
	treemat = np.delete(treemat,inds[1],axis=1)
	if limitinter >= 2:
			if lenmerged >= limitinter:
				del namelist[inds[0]]
				treemat = np.delete(treemat,inds[0],axis=0)
				treemat = np.delete(treemat,inds[0],axis=1)
	
	
	while len(namelist) > 1:
		allcb = []
		for i2 in range(inds[0]):
			allcb = allcb + [(i2,inds[0])]
		for i2 in range(inds[0]+1,len(namelist)):
			allcb = allcb + [(inds[0],i2)]
		for i2 in allcb:
			fi = namelist[i2[0]]
			ai = namelist[i2[1]]
			alli = fi+ai
			sumres = funccdpconsygy(X0[:,alli],y)
			res1 = funccdpconsygy(X0[:,fi],y)
			res2 = funccdpconsygy(X0[:,ai],y)
			thisv = sumres - res1 -res2
			treemat[i2[0],i2[1]] = thisv
		inds = np.unravel_index(np.nanargmax(treemat, axis=None), treemat.shape)
		thev = treemat[inds]
		if thev < 0:
			output = mat(output)
			np.savetxt(outputfile,output,fmt='%s',delimiter='\t',newline='\n')
			return(1)
		thep = funpvaleconsygy(X0[:, namelist[inds[0]] ], X0[:, namelist[inds[1]] ], y, thev, iternum)
		out = stroutconsygy(namelist[inds[0]],namelist[inds[1]],proteinname)
		namelist[inds[0]] = namelist[inds[0]]+namelist.pop(inds[1])
		lenmerged = len(namelist[inds[0]])
		treemat = np.delete(treemat,inds[1],axis=0)
		treemat = np.delete(treemat,inds[1],axis=1)
		if limitinter >= 2:
			if lenmerged >= limitinter:
				del namelist[inds[0]]
				treemat = np.delete(treemat,inds[0],axis=0)
				treemat = np.delete(treemat,inds[0],axis=1)
		out = [out, thev, thep]
		output.append(out)
	output = mat(output)
	np.savetxt(outputfile,output,fmt='%s',delimiter='\t',newline='\n')
	return(1)

def itptfilecon(tempresult):
	import re
	inputnames=[]
	with open(tempresult, "r") as f:
		for line in f.readlines():
			line = line.strip('\n')
			line = line.split('\t')
			thename = line[0]
			thename = re.split('[:_]', thename)
			inputnames.append(thename)
	return(inputnames)
	# [int(i) for arr in inputnames for i in arr]

def synergycon2nd(datain,dataresult,outputfile,dataanno, digitnumber=2,iternum=0):
	ii = 0
	proteinname=[]
	with open(dataanno, "r") as f:
		for line in f.readlines():
			line = line.strip('\n')  
			line = line.split("\t")
			proteinname1 = line[1]
			proteinname.append(proteinname1)
			ii = ii + 1
	X0=[]
	X1=[]
	with open(datain, "r") as f:
		for line in f.readlines():
			line = line.strip('\n')  
			line = line.split("\t")
			distype = line[-1]
			line = line[:-1]
			X0.append(line)
			X1.append(distype)
	X0=mat(X0)

	[X0rows, X0cols] = X0.shape
	for line in range(X0cols):
		thestate = np.array(X0[:,line])
		thestate = thestate.flatten()
		thestate = digitcon(thestate,digitnumber)
		X0[:,line] = np.array([thestate]).T
	X0=np.array(X0)
	y=X1

	selectcols = itptfilecon(dataresult)
	output = []
	for smallselect in selectcols:
		select = [ct in smallselect for ct in proteinname]
		if (np.sum(select) < len(smallselect)):
			stre='_'
			out = stre.join(smallselect)
			output.append([out,'NA','NA'])
		else:
			select = np.array(select)
			selectcol = proteinname
			selectcol = np.array(selectcol)
			selectcol = selectcol[select]
			oneout = synergyconselect(X0[:,select],y,selectcol,digitnumber,iternum)
			output.append(oneout)
	output = mat(output)
	np.savetxt(outputfile,output,fmt='%s',delimiter='\t',newline='\n')
	



def synergyconselect(X0,y,proteinnames, digitnumber=2,iternum=0):
	ii = 0
	proteinname={}
	for proteinname1 in proteinnames:
		proteinname[str(ii)] = proteinname1
		ii = ii + 1
	str1 = '_'
	naoutput = str1.join(proteinnames)
	naoutput = [naoutput, 'NA', 'NA']
	#
	X0=np.array(X0)
	snpnum = X0.shape[1]
	#
	namelist = list(range(snpnum))
	namelist = [[ct] for ct in namelist]
	treemat = np.random.rand(snpnum,snpnum)
	treemat[:] = np.nan
	allcb = list(itertools.combinations(list(range(len(namelist))), 2))
	for i2 in allcb:
		fi = namelist[i2[0]]
		ai = namelist[i2[1]]
		alli = fi+ai
		sumres = funccdpconsygy(X0[:,alli],y)
		res1 = funccdpconsygy(X0[:,fi],y)
		res2 = funccdpconsygy(X0[:,ai],y)
		thisv = sumres - res1 -res2
		treemat[fi,ai] = thisv
	inds = np.unravel_index(np.nanargmax(treemat, axis=None), treemat.shape)
	thev = treemat[inds]
	if thev < 0:
			return naoutput
	if len(namelist) > 2:
		namelist[inds[0]]= namelist[inds[0]]+namelist.pop(inds[1])
		treemat = np.delete(treemat,inds[1],axis=0)
		treemat = np.delete(treemat,inds[1],axis=1)
	while len(namelist) > 2:
		allcb = []
		for i2 in range(inds[0]):
			allcb = allcb + [(i2,inds[0])]
		for i2 in range(inds[0]+1,len(namelist)):
			allcb = allcb + [(inds[0],i2)]
		for i2 in allcb:
			fi = namelist[i2[0]]
			ai = namelist[i2[1]]
			alli = fi+ai
			sumres = funccdpconsygy(X0[:,alli],y)
			res1 = funccdpconsygy(X0[:,fi],y)
			res2 = funccdpconsygy(X0[:,ai],y)
			thisv = sumres - res1 -res2
			treemat[i2[0],i2[1]] = thisv
		inds = np.unravel_index(np.nanargmax(treemat, axis=None), treemat.shape)
		thev = treemat[inds]
		if thev < 0:
			return naoutput
		namelist[inds[0]] = namelist[inds[0]]+namelist.pop(inds[1])
		treemat = np.delete(treemat,inds[1],axis=0)
		treemat = np.delete(treemat,inds[1],axis=1)
	fi = namelist[0]
	ai = namelist[1]
	alli = fi+ai
	sumres = funccdpconsygy(X0[:,alli],y)
	res1 = funccdpconsygy(X0[:,fi],y)
	res2 = funccdpconsygy(X0[:,ai],y)
	outputv = sumres - res1 -res2
	#if outputv < 0:
	#	return naoutput
	outputp = funpvaleconsygy(X0[:, fi], X0[:, ai], y, thev, iternum)
	out = stroutconsygy(fi,ai,proteinname)
	out = [out, outputv, outputp]
	return(out)
