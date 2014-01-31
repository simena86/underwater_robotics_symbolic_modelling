#!/usr/bin/python
############################################
#
# script to symbolically calculate 
# dynamic and kinematics of underwater 
# manipulator-vehicle system
# by simen andresen - simena86@gmail.com
#
############################################



#~~~´\/`~~~´\/`~~~´\/`~~~´\/`~~~´\/`~~~´\/`~~~~
#
#   	>=X>						>o=>
#  _____________
# |\            \         ____
# | \____________\       /    \	       ___
# | |            |      /      \ _ ___|
# \ |            |_____/              |___
#  \|____________|				  
#
#						<OO<
#
#______________________________________________


from sympy import *
from kinematics import *

# define symbolic variables
phi,theta,psi = symbols('phi theta psi')
p,q,r = symbols('p q r')
omega_B=symbols('omega1:3')
pd,qd,rd = symbols('pd qd rd')
u,v,w = symbols('u v w')
ud,vd,wd = symbols('ud vd wd')
q = symbols('q0:7')
qd = symbols('\dot{q}_0:7')
qdd = symbols('\ddot{q}_0:7')
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

# mathematical operators
cross = symbols('\times')

# constant parameters
g = symbols('g')  	      	# gravity = 9.81 m/s^2
rho = symbols('rho')      	# density of water
nabla = symbols('nabla')	# displacement of water - volume

# Initialize parameters and variables 
omega_i=[]        	 # angular velocity of frame i as seen from inertial frame
omega_cross_i=[]     # angular velocity of frame i as seen from inertial frame -- skew symmetric representation
omegad_i=[]	         # angular velocity of frame i as seen from inertial frame
omegad_cross_i=[]	 # angular velocity of frame i as seen from inertial frame -- skew symmetric representation
v_i=[]	 	         # linear velocity of frame i as seen from inertial frame
v_ic=[]	 	         # linear velocity of center of mass of frame i as seen from inertial frame
a_i=[]	 			 # linear acceleration of center of mass of frame i as seen from inertial frame
R_i=[]				 # rotation matrix representing rotation from frame i-1 to i
r_i=[]				 # constant vector from the origin of frame i-1 to frame i, represented in frame i
r_ic=[]				 # vector from origin of frame i to center of mass of link i represented in frame i 
r_i-1c=[] 			 # vector from origin of frame i-1 to center of gravity, expressed in frame i
M_i=[]               # Mass matrix of link i
F_i=[]               # Inertial force
I_i=[]               # Moment of inertia of link i
T_i=[]               # Inertial moment
f_i=[]               # Total forces on link i
mu_i=[]              # Total moment on link i
m_i=[]               # Total mass of link i
g_i=[]				 # gravitational acceleration g_i = R_I^i*[0 0 9.81]^T on link i , i.e. the downwards acceleration expressed in the frame of link i
r_b_i=[]			 # vector from the origin of frame i-1 to center of boyancy, expressed in frame i


# Populate variables with Matrix symbols
for i in range(0,7):
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
	r_i-1c.append(MatrixSymbol('r_{'+str(i-1)+'c}{'+str(i)+'}',3,1))
	M_i.append(MatrixSymbol('M_'+str(i),3,3))
	F_i.append(MatrixSymbol('F_'+str(i),3,3))
	I_i.append(MatrixSymbol('I_'+str(i),3,3))
	T_i.append(MatrixSymbol('T_'+str(i),3,3))
	g_i.append(MatrixSymbol('g_'+str(i),3,1))
	r_b_i.append(MatrixSymbol('r_b_i'+str(i),3,1))



# inertial forces on the first body
i=0
F_i[i] = M_i[i]*(a_i[i] + omegad_cross_i[i]*r_ic[i] + omega_cross_i[i]*(omega_cross_i[i]*r_ic[i]))
T_i[i] = I_i[i]*omegad_i[i] + omega_cross_i[i]*(I_i[i]*omega_i[i])

# forward iteration for kinematics 
# and inertial forces on the generic bodies
for i in range(1,7):
	omega_i[i] = R_i[i]*(omega_i[i-1] + qd[i]*z)
	omegad_i[i] = R_i[i]*(omegad_i[i-1]+omega_cross_i[i-1]*qd[i]*z + qdd[i]*z)
	v_i[i] = R_i[i]*v_i[i-1] + omega_cross_i[i]*r_i[i-1]
	v_ic[i] = R_i[i]*v_i[i] + omega_cross_i[i]*r_ic[i-1]
	a_i[i] = R_i[i]*a_i[i-1] + omegad_cross_i[i]*r_i[i-1] + omega_cross_i[i]*(omega_cross_i[i]*r_i[i-1])
	F_i[i] = M_i[i]*(a_i[i] + omegad_cross_i[i]*r_ic[i] + omega_cross_i[i]*(omega_cross_i[i]*r_ic[i]))
	T_i[i] = I_i[i]*omegad_i[i] + omega_cross_i[i]*(I_i[i]*omega_i[i])


i=6
f_i[i] = F_i[i] - m_i[i]*g_i[i] + rho*nabla*g_i[i] 
mu_i[i] = r_i-1c[i]*cross*F_i[i] 

preview(mu_i[i])

























#i=1
# start iteration by assigning the first link of the manipulator chain based on vehicle states
# kinematics
#omega_i[i] = R_i[i]*(omega_b + qd[i]*z)
#omegad_i[i] = R_i[i]*(omegad_b+omega_cross_b*qd[i]*z + qdd[i]*z)
#v_i[i] = R_i[i]*v_b + omega_cross_i[i]*r_i[i-1]
#v_ic[i] = R_i[i]*v_b + omega_cross_i[i]*r_ic[i-1]
#a_i[i] = R_i[i]*a_b + omegad_cross_i[i]*r_i[i-1] + omega_cross_i[i]*(omega_cross_i[i]*r_i[i-1])
#inertial forces
#F_i[i] = M_i[i]*(a_i[i] + omegad_cross_i[i]*r_ic[i] + omega_cross_i[i]*(omega_cross_i[i]*r_ic[i]))
#T_i[i] = I_i[i]*omegad_i[i] + omega_cross_i[i]*(I_i[i]*omega_i[i])
