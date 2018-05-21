import time
import argparse 
import navio.mpu9250
import numpy as np
import flight_functions as fcn

IMU = navio.mpu9250.MPU9250()

IMU.initialize(1,0x06)#lowpass->20Hz
time.sleep(1)

m6a0=np.zeros((3,))
m6g0=np.zeros((3,))

for i in range(1,1000):
    ma0, mg0 = IMU.getMotion6()
    m6a0=m6a0+np.array(ma0)
    m6g0=m6g0+np.array(mg0)
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
Ts=1/250
Yaw=0
Tri=np.zeros((2,3))

c=np.array([[1,0],[0,1]])
q=np.array([[1.74E-3*Ts*Ts,0],[0,1.74E-3*Ts*Ts]])
b=np.array([[1,0],[0,1]])
r=np.array([[1*Ts*Ts,0],[0,1*Ts*Ts]])
Cgy=np.eye(3,3)*1.74E-3*Ts*Ts

while(True):
 t_estimate_Attitude0=time.time()
 m6a, m6g = IMU.getMotion6()
 m6a[0],m6a[1]=-m6a[1],-m6a[0]
 m6g[0],m6g[1]=m6g[1],m6g[0]

 m6a=np.array(m6a)
 m6g=np.array(m6g)-m6g0
 m6g[2]=-m6g[2]

 Tri=fcn.get_Trigonometrxic(x)
 J=fcn.Jacobian_forprocessvariance2(Tri)
 Jt=J.transpose()
 x,P=fcn.Kalman_filer2(x,fcn.get_angle_acc(m6a),m6g,c,b,q,r,P,Ts,Tri)
 Yaw=Yaw+(Tri[1,1]/Tri[0,0]*m6g[1]+Tri[1,0]/Tri[0,0]*m6g[2])*Ts
 phi=np.array([Yaw,x[0]-x0[0],x[1]-x0[1]])
 print(phi)
 time.sleep(Ts-time.time()+t_estimate_Attitude0)