#!/usr/bin/env python
import time
import argparse 
import sys
import navio.mpu9250
import navio.util
import numpy as np
import navio.rcinput
import navio.pwm
import rospy
import flight_functions as fcn
#from flight_cont_ver1.msg import Drone
from Robot_parameters import Drone_only


Robot=Drone_only()
IMU = navio.mpu9250.MPU9250()

#=================================Initial setup==============================================
#-------------ROS-----------------------
FRQ=250.0
rospy.init_node('test_control_ver2')
pub=rospy.Publisher('Datas',Drone,queue_size=50)
rate=rospy.Rate(FRQ)

Datas=Drone()
Datas.Sensor.header.frame_id = "/base_link"
#--------------------RC------------------
rcin = navio.rcinput.RCInput()
ref=np.zeros((3,)) #[Yaw Roll Pith]
Pitch_Max=25*np.pi/180
Roll_Max=25*np.pi/180
Yaw_Max=90*np.pi/180
RC_MIN=1104
RC_MAX=1924
Throtle_Max=1.5*np.pi
#---------------------------------------
#--------------IMU---------------------
IMU.initialize(1,0x06)#lowpass->20Hz
time.sleep(1)

m6a0=np.zeros((3,))
m6g0=np.zeros((3,))

for i in range(1,1000):
    ma0, mg0 = IMU.getMotion6()
    m6a0=m6a0+np.array(ma0)
    m6g0=m6g0+np.array(mg0)
#m6a0=m6a0/500
m6g0=m6g0/1000
m6a0=m6a0/1000
m6a0[0],m6a0[1]=m6a0[1],m6a0[0]
m6g0[0],m6g0[1]=m6g0[1],m6g0[0]

x0=np.zeros((2,))
for i in range(1,1000):
    m6a, m6g = IMU.getMotion6()
    m6a[0],m6a[1]=-m6a[1],-m6a[0]      
    m6a=np.array(m6a)
    x0=x0+fcn.get_angle_acc(m6a)
x0=x0/1000

x=x0
P=np.zeros((2,2))
Ts=1/FRQ
Yaw=0


c=np.array([[1,0],[0,1]])
q=np.array([[1.74E-3*Ts*Ts,0],[0,1.74E-3*Ts*Ts]])
b=np.array([[1,0],[0,1]])
r=np.array([[0.008*Ts*Ts,0],[0,0.008*Ts*Ts]])
Cgy=np.eye(3,3)*1.74E-3*Ts*Ts
#-------------------------------------
#-----------Controller------------------
p=np.zeros((3,))
U=np.zeros((4,))
Ixx=Robot.Ixx
Iyy=Robot.Iyy
Izz=Robot.Izz
a0=(Ixx-Iyy)/Izz
a1=(Izz-Ixx)/Iyy
a2=(Iyy-Izz)/Ixx

C10=4
C11=4
C12=6
C20=10
C21=6
C22=10
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
Mixing=Robot.Mixing
#-------------------------------------
#----------Output---------------------
PWM=np.zeros((4,))+1.000
motor_frq=400
#--------------------------------------
#=============================================================================================

def loop():
    t=time.time()
    #=========================Get RCInput and Traslate to reference=================
    t_RC=time.time()
    ref=np.array([float(rcin.read(0)),float(rcin.read(1)),float(rcin.read(3))])
    Throttle=int(rcin.read(2))
    ref[0]=fcn.map(ref[0],RC_MIN,RC_MAX,-Yaw_Max,Yaw_Max)
    ref[1]=fcn.map(ref[1],RC_MIN,RC_MAX,-Roll_Max,Roll_Max)
    ref[2]=fcn.map(ref[2],RC_MIN,RC_MAX,-Pitch_Max,Pitch_Max)
    Throttle=fcn.map(Throttle,RC_MIN,RC_MAX,0,Throtle_Max)
    t_RC=time.time()-t_RC
    #===============================================================================
    #=========================Estimate Attitude===================================
    t_estimate_Attitude0=time.time()
    m6a, m6g = IMU.getMotion6()
    m6a[0],m6a[1]=-m6a[1],-m6a[0]
    m6g[0],m6g[1]=m6g[1],m6g[0]
                    
    m6a=np.array(m6a)
    #m6a=fcn.lowpassfilter(np.array(m6a),m6a0)#-m6a0
    m6g=np.array(m6g)-m6g0
    m6g[2]=-m6g[2]
    J=fcn.Jacobian_forprocessvariance(x)
    Jt=J.transpose()
    q=np.dot(J,np.dot(Cgy,Jt))*Ts*Ts
    x,P=fcn.Kalman_filer(x,fcn.get_angle_acc(m6a),m6g,c,b,q,r,P,Ts)
    Yaw=Yaw+(np.sin(x[1])/np.cos(x[0])*m6g[1]+np.cos(x[1])/np.cos(x[0])*m6g[2])*Ts
    phi=np.array([Yaw,x[0]-x0[0],x[1]-x0[1]])
    #m6a0=m6a
    t_estimate_Attitude=time.time()-t_estimate_Attitude0
    #=============================================================================

    #=========================Controller============================================
    t_cont0=time.time()
    p=p+phi*Ts
    phi_dot=np.array([np.sin(x[1])/np.cos(x[0])*m6g[1]+np.cos(x[1])/np.cos(x[0])*m6g[2],\
                      np.cos(x[1])*m6g[1]-np.sin(x[1])*m6g[2],\
                      m6g[0]+np.sin(x[1])*np.tan(x[0])*m6g[1]+np.cos(x[1])*np.tan(x[0])*m6g[2]])
    U[0]=Throttle
    U[1]=Izz*(K10*(ref[0]-phi[0])-K20*phi_dot[0]+K30*p[0]-a0*phi_dot[1]*phi_dot[2])
    U[2]=Iyy*(K11*(ref[1]-phi[1])-K21*phi_dot[1]+K31*p[1]-a1*phi_dot[2]*phi_dot[0])
    U[3]=Ixx*(K12*(ref[2]-phi[2])-K22*phi_dot[2]+K32*p[2]-a2*phi_dot[0]*phi_dot[1])
    Thrusts=np.dot(Mixing,U)#0:FL,1:FR,2:RR,3:RL
    t_cont=time.time()-t_cont0
    #===============================================================================
    #=========================Output================================================
    t_output0=time.time()
    PWM=np.array([fcn.thrust2pwm(Thrusts[j]) for j in range(0,4)])
    #PWM[1]=fcn.map(Throttle,RC_MIN,RC_MAX,1.0/8,2.1/8)
    motor1.set_duty_cycle(PWM[0])
    motor2.set_duty_cycle(PWM[1])
    motor3.set_duty_cycle(PWM[2])
    motor4.set_duty_cycle(PWM[3])
    t_output=time.time()-t_output0
    #===============================================================================
    #=========================Debug==========================================================
    t_debag0=time.time()                    
    Datas.Sensor.linear_acceleration.x=m6a[0]
    Datas.Sensor.linear_acceleration.y=m6a[1]
    Datas.Sensor.linear_acceleration.z=m6a[2]
    Datas.Sensor.angular_velocity.x=m6g[0]
    Datas.Sensor.angular_velocity.y=m6g[1]
    Datas.Sensor.angular_velocity.z=m6g[2]
    Datas.Sensor.header.stamp=rospy.Time.now()
    Datas.Euler.x=phi[2]
    Datas.Euler.y=phi[1]
    Datas.Euler.z=phi[0]
    Datas.Inputs.FL=Thrusts[0]
    Datas.Inputs.FR=Thrusts[1]
    Datas.Inputs.RR=Thrusts[2]
    Datas.Inputs.RL=Thrusts[3]
    Datas.Ref.x=ref[2]
    Datas.Ref.y=ref[1]
    Datas.Ref.z=ref[0]
    Datas.Trottle=Throttle
    pub.publish(Datas)
    t_debug=time.time()-t_debag0
    #print "RC:", "{:7.5f}".format(t_RC),"ANG:", "{:7.5f}".format(t_estimate_Attitude),"CNT:", "{:7.5f}".format(t_cont),"OUT:", "{:7.5f}".format(t_output),
    #print "ALL", "{:7.5f}".format(time.time()-t)
    print PWM
    #=========================================================================================

    
