#!/usr/bin/env python
import time
import argparse 
import sys
import numpy as np
import rospy
import flight_functions as fcn
from std_msgs.msg import Float64MultiArray
from Robot_parameters import Drone_only
from original_msgs.msg import Droneinfo
from geometry_msgs.msg import Vector3
import tf

Robot=Drone_only()


#=================================Initial setup==============================================
#-------------ROS-----------------------
FRQ=250.0
rospy.init_node('simulation_atitude_control')
pub1=rospy.Publisher('inputs',Float64MultiArray,queue_size=50)
pub2=rospy.Publisher('passive_odom',Vector3,queue_size=50)
Inputs=Float64MultiArray()
orien=Vector3()
#-------------------------------------
#-----------Controller------------------
Ts=0.001
ref=np.zeros((3,))

p=np.zeros((3,))
U=np.zeros((4,))
Ixx=Robot.Ixx
Iyy=Robot.Iyy
Izz=Robot.Izz
a0=(Ixx-Iyy)/Izz
a1=(Izz-Ixx)/Iyy
a2=(Iyy-Izz)/Ixx

C10=600
C11=12
C12=12
C20=15
C21=3
C22=3
ram0=0.01
ram1=0.01
ram2=0.01

K10=(1+ram0+C10*C20)
K11=(1+ram1+C11*C21)
K12=(1+ram2+C12*C22)
K20=C10+C20
K21=C11+C21
K22=C12+C22
K30=ram0*C20
K31=ram1*C21
K32=ram2*C22

IMAX=0.1
IMIN=-0.1

Mixing=Robot.Mixing
t0=0
#-------------------------------------
#=============================================================================================

def callback(msg,p):

    P=np.array([msg.odom.pose.pose.position.x,msg.odom.pose.pose.position.y,msg.odom.pose.pose.position.z])
    P_dot=np.array([msg.odom.twist.twist.linear.x,msg.odom.twist.twist.linear.y,msg.odom.twist.twist.linear.z])
    phi=np.array([msg.odom.pose.pose.orientation.x,msg.odom.pose.pose.orientation.y,msg.odom.pose.pose.orientation.z,msg.odom.pose.pose.orientation.w])
    phi=tf.transformations.euler_from_quaternion(phi,'rxyz')
    phi=np.array(phi)
    phi[0],phi[2]=phi[2],phi[0]
    phi_dot=np.array([msg.odom.twist.twist.angular.x,msg.odom.twist.twist.angular.y,msg.odom.twist.twist.angular.z])
    phi_dot[0],phi_dot[2]=phi_dot[2],phi_dot[0]
    load=msg.load.force.z
    #phi[2]=-phi[2]
    phi[1]=-phi[1]
    #phi_dot[2]=-phi_dot[2]
    #phi_dot[1]=-phi_dot[1]
    p=p+phi*Ts
    Ux=0.1*(0.0-P[0])-0.02*P_dot[0]
    Uy=0.1*(0.0-P[1])-0.02*P_dot[1]

    U[0]=Robot.M*9.81#*np.cos(phi[1])
    ref[2]=np.arcsin((-Ux/U[0]*np.sin(ref[0])+Uy/U[0]*np.cos(ref[0])))
    ref[1]=-np.arcsin((-Ux/U[0]*np.cos(ref[0])-Uy/U[0]*np.sin(ref[0]))/np.cos(-ref[2]))
    U[1]=Izz*(K10*(ref[0]-phi[0])-K20*phi_dot[0]+K30*p[0]-a0*phi_dot[1]*phi_dot[2])
    U[2]=Iyy*(K11*(ref[1]-phi[1])-K21*phi_dot[1]+K31*p[1]-a1*phi_dot[2]*phi_dot[0])
    U[3]=Ixx*(K12*(ref[2]-phi[2])-K22*phi_dot[2]+K32*p[2]-a2*phi_dot[0]*phi_dot[1])
    #U[1]=-U[1]
    U[2]=-U[2]
    Thrusts=np.dot(Mixing,U)
    Inputs.data=Thrusts
    orien.x=phi[2]
    orien.y=phi[1]
    orien.z=phi[0]
    print ref[1],ref[2]
    pub1.publish(Inputs)
    pub2.publish(orien)

sub=rospy.Subscriber('info',Droneinfo,callback,p)
rospy.spin()