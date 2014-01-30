#!/usr/bin/python

from sympy import *

def RotZ(psi):
	mat = Matrix([[cos(psi), -sin(psi), 0 ],[sin(psi), cos(psi), 0 ],[0,0,1]])
	return mat

def RotX(phi):
	mat=Matrix([[1,0,0],[0, cos(phi), -sin(phi)],[ 0 , sin(phi), cos(phi)]])
	return mat 

def omegaVec(postfix):
	astr='p_' + str(postfix) + ' q_' + str(postfix) + ' r_' + str(postfix) 
	p,q,r=symbols(astr)
	mat=Matrix([[p],[q],[r]])
	return mat

def omegadVec(postfix):
	str='\dot{p}_' + postfix + ' \dot{q}_' + postfix + ' \dot{r}_' + postfix 
	p,q,r=symbols(str)
	mat=Matrix([[p],[q],[r]])
	return mat

#def vVec(postfix):
#	astr='u_' + str(postfix) + ' v_' + str(postfix) + ' w_' + str(postfix) 
#	p,q,r=symbols(astr)
#	mat=Matrix([[p],[q],[r]])

#def aVec(postfix):
#	astr='\dot{u}_' + str(postfix) + ' \dot{v}_' + str(postfix) + ' \dot{w}_' + str(postfix) 
#	p,q,r=symbols(astr)
#	mat=Matrix([[p],[q],[r]])

#def skew(vec):
#	mat= Matrix([[0 , -vec[2] ,vec[1]],[vec[2], 0, -vec[0]],[ -vec[1], vec[0], 0 ]])
#	return mat


def vec3(name,superscript,subscript):
	mat=Matrix([[var(name + '_{'+subscript+'}^{' +superscript+ '}' )],[ var(name + '_{'+subscript+'}^{' +superscript+ '}' )],[  var(name + '_{'+subscript+'}^{' +superscript+ '}' )]])
	#preview(eye(3)*mat)
	return mat

