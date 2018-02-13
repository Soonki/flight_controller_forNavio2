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
from sensor_msgs.msg import Imu
from std_msgs.msg import Float64
from geometry_msgs.msg import Vector3
import flight_functions as fcn
#from flight_cont_ver1.msg import InputForces_drone

navio.util.check_apm()

parser = argparse.ArgumentParser()
parser.add_argument("-i", help = "Sensor selection: -i [sensor name]. Sensors names: mpu is MPU9250, lsm is LSM9DS1")

if len(sys.argv) == 1:
    print "Enter parameter"
    parser.print_help()
    sys.exit(1)
elif len(sys.argv) == 2:
    sys.exit("Enter sensor name: mpu or lsm")

args = parser.parse_args()

if args.i == 'mpu':
    print "Selected: MPU9250"
    IMU = navio.mpu9250.MPU9250()
elif args.i == 'lsm':
    print "Selected: LSM9DS1"
    IMU = navio.lsm9ds1.LSM9DS1()
else:
    print "Wrong sensor name. Select: mpu or lsm"
    sys.exit(1)



if IMU.testConnection():
    print "Connection established: True"
else:
    sys.exit("Connection established: False")


#=================================Initial setup==============================================
#-------------ROS-----------------------
FRQ=250.0
rospy.init_node('test_control_ver1')
pub_imu=rospy.Publisher('imu',Imu,queue_size=50)
pub_ref_ang=rospy.Publisher('ref_ang',Vector3,queue_size=50)
pub_euler=rospy.Publisher('euler',Vector3,queue_size=50)
pub_throttle=rospy.Publisher('throttle',Float64,queue_size=50)
pub_thrusts=rospy.Publisher('thrusts',InputForces_drone,queue_size=50)
rate=rospy.Rate(FRQ)

imu=Imu()
imu.header.frame_id = "/base_link"
euler=Vector3()
thrusts=InputForces_drone()
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
IMU.initialize(1,0x04)#lowpass->20Hz
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

x=np.zeros((2,))
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
Ixx=0.0001
Iyy=0.0001
Izz=0.0002
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
Mixing=np.eye(4,4)#np.array([])
#-------------------------------------
#----------Output---------------------
PWM=np.zeros((4,))+1.000
motor_frq=50

#--------------------------------------
#=============================================================================================

with navio.pwm.PWM(0) as motor1:
    motor1.set_period(motor_frq)
    motor1.enable()
    with navio.pwm.PWM(1) as motor2:
        motor2.set_period(motor_frq)
        motor2.enable()
        with navio.pwm.PWM(2) as motor3:
            motor3.set_period(motor_frq)
            motor3.enable()
            with navio.pwm.PWM(3) as motor4:
                motor4.set_period(motor_frq)
                motor4.enable()

                while not rospy.is_shutdown():
                    t=time.time()

                    #=========================Get RCInput and Traslate to reference=================
                    t_RC=time.time()
                    ref=np.array([float(rcin.read(0)),float(rcin.read(3)),float(rcin.read(1))])
                    Throtle=int(rcin.read(2))
                    ref[0]=fcn.map(ref[0],RC_MIN,RC_MAX,-Yaw_Max,Yaw_Max)
                    ref[1]=fcn.map(ref[1],RC_MIN,RC_MAX,-Roll_Max,Roll_Max)
                    ref[2]=-fcn.map(ref[2],RC_MIN,RC_MAX,-Pitch_Max,Pitch_Max)
                    Throtle=fcn.map(Throtle,RC_MIN,RC_MAX,0,Throtle_Max)
                    t_RC=time.time()-t_RC
                    #===============================================================================

                    #=========================Estimate Attitude===================================
                    t_estimate_Attitude0=time.time()
                    m6a, m6g = IMU.getMotion6()
                    m6a=fcn.lowpassfilter(np.array(m6a),m6a0)#-m6a0
                    m6g=np.array(m6g)-m6g0
                    J=fcn.Jacobian_forprocessvariance(x)
                    Jt=J.transpose()
                    q=np.dot(J,np.dot(Cgy,Jt))*Ts*Ts
                    x,P=fcn.Kalman_filer(x,fcn.get_angle_acc(m6a),m6g,c,b,q,r,P,Ts)
                    Yaw=Yaw+(np.sin(x[1])/np.cos(x[0])*m6g[1]+np.cos(x[1])/np.cos(x[0])*m6g[2])*Ts
                    phi=np.array([Yaw,x[0],x[1]])
                    m6a0=m6a
                    t_estimate_Attitude=time.time()-t_estimate_Attitude0
                    #=============================================================================

                    #=========================Controller============================================
                    t_cont0=time.time()
                    p=p+phi*Ts
                    phi_dot=np.array([np.sin(x[1])/np.cos(x[0])*m6g[1]+np.cos(x[1])/np.cos(x[0]),np.cos(x[1])*m6g[1]-np.sin(x[1])*m6g[2],m6g[0]+np.sin(x[1])*np.tan(x[0])*m6g[1]+np.cos(x[1])*np.tan(x[0])*m6g[2]])
                    U[0]=3#Throtle
                    U[1]=Izz*(K10*(ref[0]-phi[0])-K20*phi_dot[0]+K30*p[0]-a0*phi_dot[1]*phi_dot[2])
                    U[2]=Iyy*(K11*(ref[1]-phi[1])-K21*phi_dot[1]+K31*p[1]-a1*phi_dot[2]*phi_dot[0])
                    U[3]=Ixx*(K12*(ref[2]-phi[2])-K22*phi_dot[2]+K32*p[2]-a2*phi_dot[0]*phi_dot[1])
                    Thrusts=np.dot(Mixing,U)#0:FL,1:FR,2:RR,3:RL
                    t_cont=time.time()-t_cont0
                    #===============================================================================
                    #=========================Output================================================
                    t_output0=time.time()
                    PWM=abs(Thrusts)
                    #PWM[0]=map(Throtle,RC_MIN,RC_MAX,0.0,2.0)
                    motor1.set_duty_cycle(PWM[0])
                    motor2.set_duty_cycle(PWM[1])
                    motor3.set_duty_cycle(PWM[2])
                    motor4.set_duty_cycle(PWM[3])
                    t_output=time.time()-t_output0
                    #===============================================================================
                    #=========================Debug==========================================================
                    #print "Acc:", "{:+7.3f}".format(m6a[0]), "{:+7.3f}".format(m6a[1]), "{:+7.3f}".format(m6a[2]),
                    #print " Gyr:", "{:+8.3f}".format(m6g[0]), "{:+8.3f}".format(m6g[1]), "{:+8.3f}".format(m6g[2]),
                    #print " Mag:", "{:+7.3f}".format(m9m[0]), "{:+7.3f}".format(m9m[1]), "{:+7.3f}".format(m9m[2]),
                    #print " Ang:", "{:+7.3f}".format(phi[0]), "{:+7.3f}".format(phi[1]),"{:+7.3f}".format(phi[2]),
                    #print " Ref:", "{:+7.3f}".format(ref[0]), "{:+7.3f}".format(ref[1]), "{:+7.3f}".format(ref[2]),
                    #print " U:", "{:+7.3f}".format(U[0]), "{:+7.3f}".format(U[1]), "{:+7.3f}".format(U[2]), "{:+7.3f}".format(U[3]),
                    t_debag0=time.time()
                    
                    imu.linear_acceleration.x=m6a[0]
                    imu.linear_acceleration.y=m6a[1]
                    imu.linear_acceleration.z=m6a[2]
                    imu.angular_velocity.x=m6g[0]
                    imu.angular_velocity.y=m6g[1]
                    imu.angular_velocity.z=m6g[2]
                    quataniom=fcn.EulertoQuaternion(phi[2],phi[1],phi[0])
                    imu.orientation.w=quataniom[0]
                    imu.orientation.z=quataniom[1]
                    imu.orientation.y=quataniom[2]
                    imu.orientation.x=quataniom[3]
                    imu.header.stamp=rospy.Time.now()
                    euler.x=phi[2]
                    euler.y=phi[1]
                    euler.z=phi[0]
                    thrusts.FL=Thrusts[0]
                    thrusts.FR=Thrusts[1]
                    thrusts.RR=Thrusts[2]
                    thrusts.RL=Thrusts[3]
                    pub_euler.publish(euler)
                    pub_imu.publish(imu)
                    #pub_ref_ang.publish(ref)
                    pub_throttle.publish(Throtle)
                    pub_thrusts.publish(thrusts)
                    t_debug=time.time()-t_debag0
                    print "RC:", "{:7.5f}".format(t_RC),"ANG:", "{:7.5f}".format(t_estimate_Attitude),"CNT:", "{:7.5f}".format(t_cont),"OUT:", "{:7.5f}".format(t_output),
                    print time.time()-t
                    #=========================================================================================
                    rate.sleep()
