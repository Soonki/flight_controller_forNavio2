#!/usr/bin/env python

import numpy as np

#================================================Function of Kalmanfilter====================================================
#state-equations
#x=[ph th ps]
#xk+1=f(xk)+bw
#yk=cxk+v
#w->q v->r
def lowpassfilter(x,y):
    return 0.001*y+0.999*x

def get_angle_acc(acc):
    th_acc=np.arctan2(-acc[0],np.sqrt(acc[1]*acc[1]+acc[2]*acc[2]))
    ps_acc=np.arctan2(acc[1],acc[2])
    y=np.array([th_acc,ps_acc])
    return y

def get_preEstimation(x,gyro,Ts):
    Q=np.array([[0,np.cos(x[1]),-np.sin(x[1])],[1,np.sin(x[1])*np.tan(x[0]),np.cos(x[1])*np.tan(x[0])]])
    return x+np.dot(Q,gyro)*Ts

def get_preVariance(x,gyro,P,b,q,Ts):
    A=np.array([[1,-(np.sin(x[1])*gyro[1]+np.cos(x[1])*gyro[2])*Ts],[(np.cos(x[1])/np.cos(x[0])/np.cos(x[0])*gyro[1]-np.sin(x[1])/np.cos(x[0])/np.cos(x[0])*gyro[2])*Ts,1+(np.cos(x[1])*np.tan(x[0])*gyro[1]-np.sin(x[1])*np.tan(x[0])*gyro[2])*Ts]])
    At=A.transpose()
    return np.dot(A,np.dot(P,At))+q

def get_Kalamgain(P,c,r):
    CPC=np.dot(c,np.dot(P,c))+r
    return np.dot(P,np.dot(c,np.linalg.inv(CPC)))

def get_Variance(g,c,P):
    return np.dot(np.array([[1,0],[1,0]])-np.dot(g,c),P)

def Kalman_filer(x,y,gyro,c,b,q,r,P,Ts):
    x_=get_preEstimation(x,gyro,Ts)
    P_=get_preVariance(x,gyro,P,b,q,Ts)
    g=get_Kalamgain(P_,c,r)
    return x_+np.dot(g,y-np.dot(c,x_)),get_Variance(g,c,P_)

def Jacobian_forprocessvariance(x):
    Cps=np.cos(x[1])
    Sps=np.sin(x[1])
    Tth=np.tan(x[0])
    return np.array([[0,Cps,-Sps],[1,Sps*Tth,Cps*Tth]])

#------------------------------------------------------------------------------

def get_preEstimation2(x,gyro,Ts,Tri):
    Q=np.array([[0,Tri[1,0],-Tri[1,1]],[1,Tri[1,1]*Tri[0,2],Tri[1,0]*Tri[0,2]]])
    return x+np.dot(Q,gyro)*Ts

def get_preVariance2(x,gyro,P,b,q,Ts,Tri):
    A=np.array([[1,-(Tri[1,1]*gyro[1]+Tri[1,0]*gyro[2])*Ts],[(Tri[1,0]/Tri[0,0]/Tri[0,0]*gyro[1]-Tri[1,1]/Tri[0,0]/Tri[0,0]*gyro[2])*Ts,1+(Tri[1,0]*Tri[0,2]*gyro[1]-Tri[1,1]*Tri[0,2]*gyro[2])*Ts]])
    return np.dot(A,np.dot(P,A.transpose()))+q

def Kalman_filer2(x,y,gyro,c,b,q,r,P,Ts,Tri):
    x_=get_preEstimation2(x,gyro,Ts,Tri)
    P_=get_preVariance2(x,gyro,P,b,q,Ts,Tri)
    g=get_Kalamgain(P_,c,r)
    return x_+np.dot(g,y-np.dot(c,x_)),get_Variance(g,c,P_)

def Jacobian_forprocessvariance2(Tri):
    return np.array([[0,Tri[1,0],-Tri[1,1]],[1,Tri[1,1]*Tri[0,2],Tri[1,0]*Tri[0,2]]])

#===============================================================================================

#======================================Other Functions==========================================
def suturation(value,Max,Min):
    if value>Max:
        value=Max
    elif value<Min:
        value=Min

def map(x,x_min,x_max,Min,Max):
    return (Max-Min)/(x_max-x_min)*(x-x_min)+Min

def thrust2pwm(T):
    if T<=0:
        return 0.1
    a=209.8
    b=-17.84
    c=-2.004
    pwm=(-b+np.sqrt(b*b-4*a*(c-T)))/2/a
    return pwm

def get_Trigonometrxic(x):
    return np.array([[np.cos(x[0]),np.sin(x[0]),np.tan(x[0])],[np.cos(x[1]),np.sin(x[1]),np.tan(x[1])]])

def EulertoQuaternion(ph,th,ps):
    Cph=np.cos(ph*0.5)
    Sph=np.sin(ph*0.5)
    Cth=np.cos(th*0.5)
    Sth=np.sin(th*0.5)
    Cps=np.cos(ps*0.5)
    Sps=np.sin(ps*0.5)
    w=Cps*Cth*Cph+Sps*Sth*Sph
    x=Sps*Cth*Cph-Cps*Sth*Sph
    y=Cps*Sth*Cph+Sps*Cth*Sph
    z=Cps*Cth*Sph-Sps*Sth*Cph
    return np.array([w,x,y,z])
#==============================================================================================