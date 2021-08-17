# -*- coding: utf-8 -*-
"""
Created on Mon May 31 01:29:05 2021

@author: zhang
"""

import sys
import numpy as np
from numpy import mat
import itertools
import pandas as pd

def itptfileppi(tempresult):
	import re
	inputnames=[]
	with open(tempresult, "r") as f:
		for line in f.readlines():
			line = line.strip('\n')
			thename = line[0]
			thename = re.split('[:_]', thename)
			inputnames.append(thename)
	return(inputnames)


def getPPIdatain(inputnames,refname):
	reframe=[]
	with open(refname, "r") as f:
		for line in f.readlines():
			line = line.strip('\n')
			line = line.split("_")
			reframe = reframe + line
	reframe=np.array(reframe)
	select = [ct in inputnames for ct in reframe]
	select = np.array(select) 
	inputnames = reframe[select]
	output = []
	for input1name in inputnames:
		temp1col = reframe == input1name
		temp1col = temp1col.astype(np.int)
		output.append(temp1col)
	output = np.mat(output)
	output = output.T
	output = (output, inputnames)
	return(output)


			


def sygyPPI(datain,refmat,geneexp,inter=1000,changethreshold=1e-19,alpha=0.7):
	refmat = np.mat(refmat)
	geneexp = np.mat(geneexp)
	datain = np.mat(datain)
	tosave = np.array((refmat.sum(0) > 0)).flatten().tolist()
	refmat  =  refmat[tosave,:][:,tosave]
	geneexp =  geneexp[tosave,:]
	datain  =  datain[tosave,:]
	##
	refmat = np.multiply(np.multiply(geneexp,refmat),geneexp.T)
	norA = 1/ np.power(refmat.sum(1),0.5)
	refmat = np.multiply(np.multiply(norA,refmat),norA.T)
	datain = datain.astype('float64')
	for xi in range(datain.shape[1]):
		F0 = datain[:,xi]
		Ft = F0
		for n in range(inter):
			Ft2 = alpha*np.dot(refmat,Ft) + (1-alpha)*F0
			change = np.sum(abs(Ft2 - Ft))
			if change < changethreshold:
				break
			Ft = Ft2
		datain[:,xi] = Ft2
	output = (datain)
	return(output)

def funGC(X0,fi,ai,threshold):
	datain01 = X0[:,fi]
	datain02 = X0[:,ai]
	datain = datain01 + datain02
	output = np.sum(datain >= threshold) - np.sum(datain01 >= threshold) - np.sum(datain02 >= threshold)
	return(output)


def funGCpv(X0,fi,ai,iternum,threshold):
	datain01 = X0[:,fi]
	datain02 = X0[:,ai]
	datain = datain01 + datain02
	threshold2 = np.sum(datain >= threshold)
	output = np.sum(datain >= threshold) - np.sum(datain01 >= threshold) - np.sum(datain02 >= threshold)
	if output >= 1:
		count = 0
		for n in range(iternum):
			datain1 = np.random.shuffle(datain01)
			datain = datain01 + datain02
			all = np.sum(datain >= threshold)
			if all >= threshold2:
				count = count + 1
		ratio = (count + 1) / (iternum+1)
		return(ratio)
	else:
		return(1)

def stroutppi(x,y,inputnames):
	if len(inputnames):
		x = [str(inputnames[numeric_string]) for numeric_string in x]
		y = [str(inputnames[numeric_string]) for numeric_string in y]
		str1 = '_'
		x = str1.join(x)
		y = str1.join(y)
		out = x + ':' + y
		return(out)
	else:
		str1 = '_'
		x = str1.join(x)
		y = str1.join(y)
		out = x + ':' + y
		return(out)

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


def digitselectcols(selectcols, genename):
	ii = 0
	temphash={}
	for line in genename:
		temphash[line] = ii
		ii = ii + 1
	newselectcols = []
	for line in selectcols:
		temparray = []
		for x in line:
			if x in temphash:
				temparray.append(temphash[x])
			else:
				temparray = []
				break
		newselectcols.append(temparray) 
	return(newselectcols)


	
def GCjudge(expmat,dataanno,geneexp,outputfile,dataresult,repnum=100,threshold=0.001,iternum=10, changethreshold=1e-19,alpha=0.7):
	genename = np.array(pd.read_csv(dataanno, delimiter='\t', header=None))
	genename = genename.flatten().tolist()
	selectcols = itptfilecon(dataresult)
	flattennames =  [i for arr in selectcols for i in arr]
	ppiinit = getPPIdatain(flattennames,dataanno)
	genename = ppiinit[1]
	ppiinit	= ppiinit[0]
	geneexp = np.mat(pd.read_csv(geneexp, delimiter='\t', header=None))
	expmat = np.mat(pd.read_csv(expmat, delimiter='\t', header=None))
	ppiinit = sygyPPI(ppiinit,expmat,geneexp,repnum,changethreshold,alpha)
	numberselectcols = digitselectcols(selectcols, genename)
	output = []
	for smallselect, smallnames in zip(numberselectcols, selectcols):
		if (len(smallselect) == 0):
			stre='_'
			out = stre.join(smallnames)
			output.append([out,'NA','NA'])
		else:
			oneout = GCjudgeback(ppiinit[:,smallselect],threshold, iternum,smallnames)
			output.append(oneout)
	output = mat(output)
	np.savetxt(outputfile,output,fmt='%s',delimiter='\t',newline='\n')

	

def GCjudgeback(matin, threshold = 0.001, iternum = 10, inputnames=[]):
	matin=np.mat(matin)
	matinrows = matin.shape[1]
	#
	namelist = list(range(matinrows))
	namelist = [[ct] for ct in namelist]
	treemat = np.random.rand(matinrows,matinrows)
	treemat[:] = np.nan
	allcb = list(itertools.combinations(list(range(len(namelist))), 2))
	for i2 in allcb:
		fi = namelist[i2[0]]
		ai = namelist[i2[1]]
		thisv = funGC(matin,fi,ai,threshold)
		treemat[fi,ai] = thisv
	inds = np.unravel_index(np.nanargmax(treemat, axis=None), treemat.shape)
	thev = treemat[inds]
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
			treemat[i2[0],i2[1]] = funGC(matin,fi,ai,threshold)
		inds = np.unravel_index(np.nanargmax(treemat, axis=None), treemat.shape)
		thev = treemat[inds]
		out = stroutppi(namelist[inds[0]], namelist[inds[1]], inputnames)
		namelist[inds[0]] = namelist[inds[0]]+namelist.pop(inds[1])
		treemat = np.delete(treemat,inds[1],axis=0)
		treemat = np.delete(treemat,inds[1],axis=1)
	fi = namelist[0]
	ai = namelist[1]
	outputv = funGC(matin,fi,ai,threshold)
	outputp = funGCpv(matin,fi,ai,iternum,threshold)
	out = stroutppi(fi,ai,inputnames)
	out = [out, outputv, outputp]
	return(out)

def GCana(expmat,dataanno,geneexp,outputfile,dataselect='NULLdata',repnum=100,threshold=0.001,iternum=10,limitinter=5, changethreshold=1e-19,alpha=0.7):
	if dataselect == 'NULLdata':
		genename = np.array(pd.read_csv(dataanno, delimiter='\t', header=None))
		ppiinit  = np.eye(genename.shape[0], dtype=int)
		genename = genename.flatten().tolist()
		geneexp = np.mat(pd.read_csv(geneexp, delimiter='\t', header=None))
		expmat = np.mat(pd.read_csv(expmat, delimiter='\t', header=None))
		ppiinit = sygyPPI(ppiinit,expmat,geneexp,repnum,changethreshold,alpha)
		GCanaback(ppiinit, outputfile, threshold, iternum, genename,limitinter)
	else:
		dataselect = np.array(pd.read_csv(dataselect, delimiter='\t', header=None))[:,0].tolist()
		ppiinit = getPPIdatain(dataselect,dataanno)
		genename = ppiinit[1]
		ppiinit	= ppiinit[0]
		geneexp = np.mat(pd.read_csv(geneexp, delimiter='\t', header=None))
		expmat = np.mat(pd.read_csv(expmat, delimiter='\t', header=None))
		ppiinit = sygyPPI(ppiinit,expmat,geneexp,repnum,changethreshold,alpha)
		GCanaback(ppiinit, outputfile, threshold, iternum, genename,limitinter)

def GCanaback(matin, outputfile, threshold = 0.001, iternum = 10, inputnames=[],limitinter=5):
	matin=np.mat(matin)
	matinrows = matin.shape[1]
	naoutput = (0,1)
	#
	namelist = list(range(matinrows))
	namelist = [[ct] for ct in namelist]
	treemat = np.random.rand(matinrows,matinrows)
	treemat[:] = np.nan
	allcb = list(itertools.combinations(list(range(len(namelist))), 2))
	output = []
	for i2 in allcb:
		fi = namelist[i2[0]]
		ai = namelist[i2[1]]
		thisv = funGC(matin,fi,ai,threshold)
		treemat[fi,ai] = thisv
	inds = np.unravel_index(np.nanargmax(treemat, axis=None), treemat.shape)
	thev = treemat[inds]
	thep = funGCpv(matin, namelist[inds[0]], namelist[inds[1]], iternum, threshold)
	out = stroutppi(namelist[inds[0]], namelist[inds[1]], inputnames)
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
			treemat[i2[0],i2[1]] = funGC(matin,fi,ai,threshold)
		inds = np.unravel_index(np.nanargmax(treemat, axis=None), treemat.shape)
		thev = treemat[inds]
		out = stroutppi(namelist[inds[0]], namelist[inds[1]], inputnames)
		thep = funGCpv(matin, namelist[inds[0]], namelist[inds[1]], iternum, threshold)
		out = [out, thev, thep]
		output.append(out)
		namelist[inds[0]] = namelist[inds[0]]+namelist.pop(inds[1])
		lenmerged = len(namelist[inds[0]])
		treemat = np.delete(treemat,inds[1],axis=0)
		treemat = np.delete(treemat,inds[1],axis=1)
		if limitinter >= 2:
			if lenmerged >= limitinter:
				del namelist[inds[0]]
				treemat = np.delete(treemat,inds[0],axis=0)
				treemat = np.delete(treemat,inds[0],axis=1)
	output = mat(output)
	np.savetxt(outputfile,output,fmt='%s',delimiter='\t',newline='\n')


