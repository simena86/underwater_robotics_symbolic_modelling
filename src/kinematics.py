#!/usr/bin/python

from sympy import *

def RotZ(psi):
	return Matrix([[cos(psi), -sin(psi), 0 ],[sin(psi), cos(psi), 0 ],[0,0,1]])

def RotX(phi):
	return Matrix([[1,0,0],[0, cos(phi), -sin(phi)],[ 0 , sin(phi), cos(phi)]])
