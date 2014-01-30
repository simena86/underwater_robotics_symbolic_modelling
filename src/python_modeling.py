#!/usr/bin/python
############################################
#
# script to symbolically calculate 
# dynamic and kinematics of underwater 
# manipulator-vehicle system
#
############################################

from sympy import *
from kinematics import *

# define symbolic variables
phi,theta,psi = symbols('phi theta psi')
p,q,r = symbols('p q r')
omega_B=symbols('omega1:3')
pd,qd,rd = symbols('pd qd rd')
u,v,w = symbols('u v w')
ud,vd,wd = symbols('ud vd wd')
q = symbols('q0:6')
qd = symbols('\dot{q}_0:6')
qdd = symbols('\ddot{q}_0:6')

print qdd

v_b= MatrixSymbol('v_b',3,1)
v_cb= MatrixSymbol('v_{c,b}',3,1)
a_b= MatrixSymbol('a_b',3,1)


omega_b = MatrixSymbol('\omega_b',3,1)
omegad_b = MatrixSymbol('\dot{\omega}_b',3,1)
omega_cross_b = MatrixSymbol('\widehat{\omega}_b',3,3)
omegad_cross_b = MatrixSymbol('\widehat{\dot{\omega}}_b',3,3)


# DH parameters
a_DH=symbols('a_DH1:6')
alpha_DH=symbols('alpha1:6')
theta_DH=symbols('theta_DH1:6')
d_DH=symbols('d_DH1:6')

# frames
z=MatrixSymbol('z_i',3,1)

# number of links
n=6		

#initialize vectors to hold the motion of the 
omega_i=[]   # angular velocity of frame i as seen from inertial frame
omega_cross_i=[]   # angular velocity of frame i as seen from inertial frame
omegad_i=[]	 # angular velocity of frame i as seen from inertial frame
omegad_cross_i=[]	 # angular velocity of frame i as seen from inertial frame
v_i=[]	 	 # linear velocity of frame i as seen from inertial frame
v_ic=[]	 	 # linear velocity of center of mass of frame i as seen from inertial frame
a_i=[]	 	 # linear acceleration of center of mass of frame i as seen from inertial frame
R_i=[]	 	 
r_i=[]
r_ic=[]
M_i=[]
F_i=[]
I_i=[]
T_i=[]


for i in range(0,6):
	omega_i.append(MatrixSymbol('\omega_'+str(i),3,1))
	omega_cross_i.append(MatrixSymbol('\widehat{\omega}_'+str(i),3,3))
	omegad_i.append(MatrixSymbol('\dot{\omega}_'+str(i),3,1))
	omegad_cross_i.append(MatrixSymbol('\widehat{\dot{\omega}}_'+str(i),3,3))
	v_i.append(MatrixSymbol('v_'+str(i),3,1))
	v_ic.append(MatrixSymbol('v_{'+str(i)+',c}',3,1))
	a_i.append(MatrixSymbol('a_'+str(i),3,1))
	R_i.append(MatrixSymbol('R_'+str(i),3,3))
	r_i.append(MatrixSymbol('r_'+str(i),3,1))
	r_ic.append(MatrixSymbol('r_{'+str(i)+'c}',3,1))

	M_i.append(MatrixSymbol('M_'+str(i),3,3))
	F_i.append(MatrixSymbol('F_'+str(i),3,3))
	I_i.append(MatrixSymbol('I_'+str(i),3,3))
	T_i.append(MatrixSymbol('T_'+str(i),3,3))






# start iteration by assigning the first link of the manipulator chain based on vehicle states
i=1
# kinematics
omega_i[i] = R_i[i]*(omega_b + qd[i]*z)
omegad_i[i] = R_i[i]*(omegad_b+omega_cross_b*qd[i]*z + qdd[i]*z)
v_i[i] = R_i[i]*v_b + omega_cross_i[i]*r_i[i-1]
v_ic[i] = R_i[i]*v_b + omega_cross_i[i]*r_ic[i-1]
a_i[i] = R_i[i]*a_b + omegad_cross_i[i]*r_i[i-1] + omega_cross_i[i]*(omega_cross_i[i]*r_i[i-1])
#inertial forces
F_i[i] = M_i[i]*(a_i[i] + omegad_cross_i[i]*r_ic[i] + omega_cross_i[i]*(omega_cross_i[i]*r_ic[i]))
T_i[i] = I_i[i]*omegad_i[i] + omega_cross_i[i]*(I_i[i]*omega_i[i])



