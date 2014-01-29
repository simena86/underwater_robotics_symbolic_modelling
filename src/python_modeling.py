#!/usr/bin/python
#
# script to symbolically calculate 
# dynamic and kinematics of underwater 
# manipulator-vehicle system

from sympy import *
from kinematics import *

# define symbolic variables
phi,theta,psi = symbols('phi theta psi')
p,q,r = symbols('p q r')
pd,qd,rd = symbols('pd qd rd')
u,v,w = symbols('u v w')
ud,vd,wd = symbols('ud vd wd')
q = symbols('q1:6')
qd = symbols('qd1:6')
qdd = symbols('qdd1:6')
omega=symbols('omega1:6')

# DH parameters
a_DH=symbols('a_DH1:6')
alpha_DH=symbols('alpha1:6')
theta_DH=symbols('theta_DH1:6')
d_DH=symbols('d_DH1:6')

n=6		# number of links

omega[1] = RotZ(q[1])*

