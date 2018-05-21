#!/usr/bin/env python
import numpy as np
import tf
import matlab.engine

class Drone_only():
    def __init__(self):
        self.M=0.612#[kg]
        self.L=0.29#[m]
        self.l=0.275#length of the link
        self.lg=0.08
        self.m_link=0.001953
        self.Lfh=0.0905
        self.Lrh=0.0985
        self.Lfv=0.1135
        self.Lrv=0.1055
        self.Ixx=0.002496#[kgm^2]
        self.Iyy=0.00252
        self.Izz=0.004265
        self.C=10
        self.Mixing=np.array([[      1,     1,      1,     1],\
                              [self.C,-self.C,self.C,-self.C],\
                              [self.Lfh,self.Lfh,-self.Lrh,-self.Lrh],\
                              [self.Lfv,-self.Lfv,-self.Lrv,self.Lrv]])
        self.Mixing2=np.array([[      1,     1,      1,     1],\
                              [self.C,-self.C,self.C,-self.C],\
                              [self.Lfh,self.Lfh,-self.Lfh,-self.Lfh],\
                              [self.Lfh,-self.Lfh,-self.Lfh,self.Lfh]])
        self.Mixing=np.linalg.inv(self.Mixing)
        self.Mixing2=np.linalg.inv(self.Mixing2)


    def euler_from_quaternion(self,q):
        ps=np.arctan2(2.0*(q[1]*q[2]+q[3]*q[0]),1.0-2.0*(q[0]*q[0]+q[1]*q[1]))
        th=np.arcsin(-2.0*(q[0]*q[2]-q[3]*q[1]))
        ph=np.arctan2(2.0*(q[0]*q[1]+q[3]*q[2]),1.0-2.0*(q[1]*q[1]+q[2]*q[2]))
        return np.array([ph,th,ps])

    def updateSinCosTan(self):
        self.Cph,self.Cth,self.Cps=np.cos(self.phi[0]),np.cos(self.phi[1]),np.cos(self.phi[2])
        self.Sph,self.Sth,self.Sps=np.sin(self.phi[0]),np.sin(self.phi[1]),np.sin(self.phi[2])
        self.Tph,self.Tth,self.Tps=np.tan(self.phi[0]),np.tan(self.phi[1]),np.tan(self.phi[2])

    def updateQuaternion(self,msg):
        self.quaternion=np.array([msg.odom.pose.pose.orientation.x,msg.odom.pose.pose.orientation.y,msg.odom.pose.pose.orientation.z,msg.odom.pose.pose.orientation.w])
        return self.quaternion

    def updatePhi(self,msg):

        self.phi=np.array(self.euler_from_quaternion(self.updateQuaternion(msg)))
        #self.phi[0],self.phi[2]=self.phi[2],self.phi[0]

    def updatePosition(self,msg):
        self.pos=np.array([msg.odom.pose.pose.position.x,msg.odom.pose.pose.position.y,msg.odom.pose.pose.position.z])
    
    def updateVelocity(self,msg):
        self.pos_dot=np.array([msg.odom.twist.twist.linear.x,msg.odom.twist.twist.linear.y,msg.odom.twist.twist.linear.z])

    def updateAngularvelocity(self,msg):
        self.omega=np.array([msg.odom.twist.twist.angular.x,msg.odom.twist.twist.angular.y,msg.odom.twist.twist.angular.z])
    
    def getQuaternion_matrix(self,q):
        return tf.transformations.quaternion_matrix(q)
    
    def updateExternalforce(self,msg):
        self.R=self.getQuaternion_matrix(self.updateQuaternion(msg))
        R=self.R.T
        self.force=np.array([msg.load.force.x,msg.load.force.y,msg.load.force.z,msg.load.torque.x,msg.load.torque.y,msg.load.torque.z])
        fg=np.dot(R[0:3,0:3],np.array([0,0,-9.81*self.m_link]))
        tg=np.array([-self.lg*fg[1],self.lg*fg[0],0])
        self.force[0:3]=-(self.force[0:3]-fg)
        self.force[3:6]=-(self.force[3:6]-tg)

    def updateExternalforce_body(self,msg):
        self.R=self.getQuaternion_matrix(self.updateQuaternion(msg))
        R=self.R.T
        self.force=np.array([msg.load.force.x,msg.load.force.y,msg.load.force.z,msg.load.torque.x,msg.load.torque.y,msg.load.torque.z])
        fg=np.dot(R[0:3,0:3],np.array([0,0,-9.81*self.m_link]))
        tg=np.array([-self.lg*fg[1],self.lg*fg[0],0])
        self.force[0:3]=-(self.force[0:3]-fg)
        self.force[3:6]=-(self.force[3:6]-tg)

    
    def updateQmatrix(self):
        #phi_dot=Q*Angular_velocity
        self.Q=np.array([[0,self.Sps/self.Cth,self.Cps/self.Cth],\
                         [0,self.Cps,-self.Sps],\
                         [1,self.Sps*self.Tth,self.Cps*self.Tth]])

    def updatePhi_dot(self,msg):
        self.updateAngularvelocity(msg)
        self.updateQmatrix()
        self.phi_dot=np.dot(self.Q,self.omega)

    def update(self,msg):
        self.updatePhi(msg)
        self.updateSinCosTan()
        self.updatePosition(msg)
        self.updateVelocity(msg)
        self.updatePhi_dot(msg)
        self.updateExternalforce(msg)

#=============================================================================================

class Drone_one_joint():
    def __init__(self):
        self.m=0.612#[kg]
#        self.L=0.29#[m]
#        self.l=0.275#length of the link
#        self.m_link0=0.005
#        self.m_link1=0.016
#        self.m_ring=0.001953
#        self.M=self.m+self.m_link0+self.m_link1+self.m_ring
#        self.l_link0g=0.0325
#        self.l_link1g=0.086
#        self.l_ringg=0.045

#        self.l_link0=0.065
#        self.l_link1=0.086*2
#        self.l_ring=0.045

        self.m_link0=0.007019
        self.m_link1=0.021058
        self.m_ring=0.001343
        self.M=self.m+self.m_link0+self.m_link1+self.m_ring
        self.l_link0g=0.025
        self.l_link1g=0.075
        self.l_ringg=0.045

        self.l_link0=0.05
        self.l_link1=0.15
        self.l_ring=0.045

        self.Lfh=0.0905
        self.Lrh=0.0985
        self.Lfv=0.1135
        self.Lrv=0.1055
        self.Ixx=0.002496#[kgm^2]
        self.Iyy=0.00252
        self.Izz=0.004265
        self.C=10
        self.Mixing=np.array([[      1,     1,      1,     1],\
                              [self.C,-self.C,self.C,-self.C],\
                              [self.Lfh,self.Lfh,-self.Lrh,-self.Lrh],\
                              [self.Lfv,-self.Lfv,-self.Lrv,self.Lrv]])
        self.Mixing2=np.array([[      1,     1,      1,     1],\
                              [self.C,-self.C,self.C,-self.C],\
                              [self.Lfh,self.Lfh,-self.Lfh,-self.Lfh],\
                              [self.Lfh,-self.Lfh,-self.Lfh,self.Lfh]])
        self.Mixing=np.linalg.inv(self.Mixing)
        self.Mixing2=np.linalg.inv(self.Mixing2)

        self.pos=np.zeros((3,))
        self.pos_dot=np.zeros((3,))
        self.phi=np.zeros((3,))
        self.phi_dot=np.zeros((3,))
        self.joint=0.0
        self.joint_dot=0.0
        self.force=np.zeros((6,))

    def euler_from_quaternion(self,q):
        ps=np.arctan2(2.0*(q[1]*q[2]+q[3]*q[0]),1.0-2.0*(q[0]*q[0]+q[1]*q[1]))
        th=np.arcsin(-2.0*(q[0]*q[2]-q[3]*q[1]))
        ph=np.arctan2(2.0*(q[0]*q[1]+q[3]*q[2]),1.0-2.0*(q[1]*q[1]+q[2]*q[2]))
        return np.array([ph,th,ps])

    def updateSinCosTan(self):
        self.Cph,self.Cth,self.Cps=np.cos(self.phi[0]),np.cos(self.phi[1]),np.cos(self.phi[2])
        self.Sph,self.Sth,self.Sps=np.sin(self.phi[0]),np.sin(self.phi[1]),np.sin(self.phi[2])
        self.Tph,self.Tth,self.Tps=np.tan(self.phi[0]),np.tan(self.phi[1]),np.tan(self.phi[2])

    def updateQuaternion(self,msg):
        self.quaternion=np.array([msg.odom.pose.pose.orientation.x,msg.odom.pose.pose.orientation.y,msg.odom.pose.pose.orientation.z,msg.odom.pose.pose.orientation.w])
        return self.quaternion

    def updatePhi(self,msg):

        self.phi=np.array(self.euler_from_quaternion(self.updateQuaternion(msg)))
        #self.phi[0],self.phi[2]=self.phi[2],self.phi[0]

    def updatePosition(self,msg):
        self.pos=np.array([msg.odom.pose.pose.position.x,msg.odom.pose.pose.position.y,msg.odom.pose.pose.position.z])
    
    def updateVelocity(self,msg):
        self.pos_dot=np.array([msg.odom.twist.twist.linear.x,msg.odom.twist.twist.linear.y,msg.odom.twist.twist.linear.z])

    def updateAngularvelocity(self,msg):
        self.omega=np.array([msg.odom.twist.twist.angular.x,msg.odom.twist.twist.angular.y,msg.odom.twist.twist.angular.z])
        
    def updateJointangle(self,msg):
        self.joint=msg.joint.position[0]
        self.joint_dot=msg.joint.velocity[0]
    
    def getQuaternion_matrix(self,q):
        return tf.transformations.quaternion_matrix(q)
    
    def updateExternalforce(self,msg):
        self.R=self.getQuaternion_matrix(self.updateQuaternion(msg))
        R=self.R.T
        self.force=np.array([msg.load.force.x,msg.load.force.y,msg.load.force.z,msg.load.torque.x,msg.load.torque.y,msg.load.torque.z])
        
        fg=np.dot(R[0:3,0:3],np.array([0,0,-9.81*(self.m_link0+self.m_link1+self.m_ring)]))
        #tg=np.array([-self.lg*fg[1],self.lg*fg[0],0])
        self.force[0:3]=self.force[0:3]-fg
        #self.force[1],self.force[2]=-self.force[1],-self.force[2]
        self.force[0:3]=-np.dot(self.R[0:3,0:3],self.force[0:3])
        self.R_eta=np.eye(3)
        self.R_eta[1:3,1:3]=np.array([[np.cos(self.joint),np.sin(self.joint)],[-np.sin(self.joint),np.cos(self.joint)]])
        Rd=np.dot(self.R[0:3,0:3],self.R_eta)
        self.r=self.R[0:3,2]*self.l_link0+Rd[:,2]*(self.l_link1+self.l_ring)
        self.tau=np.cross(self.r,self.force[0:3])
        #print(self.r)
        #self.force[3:6]=-(self.force[3:6]-tg)
        #print(self.force[0:3])

    
    def updateQmatrix(self):
        #phi_dot=Q*Angular_velocity
        self.Q=np.array([[0,self.Sps/self.Cth,self.Cps/self.Cth],\
                         [0,self.Cps,-self.Sps],\
                         [1,self.Sps*self.Tth,self.Cps*self.Tth]])

    def updatePhi_dot(self,msg):
        self.updateAngularvelocity(msg)
        self.updateQmatrix()
        self.phi_dot=np.dot(self.Q,self.omega)

    def update(self,msg):
        self.updatePhi(msg)
        self.updateSinCosTan()
        self.updatePosition(msg)
        self.updateVelocity(msg)
        self.updatePhi_dot(msg)
        self.updateJointangle(msg)
        self.updateExternalforce(msg)

if __name__ == '__main__':
    R=Drone_one_joint()
    q=np.array([0.3,0,0,0.2])
    q=tf.transformations.quaternion_from_euler(0.8,0.8,0.1,'rzyx')
    print tf.transformations.quaternion_matrix(q)
    print type(tf.transformations.quaternion_matrix(q))
    #print R.euler_from_quaternion(q),tf.transformations.euler_from_quaternion(q,'rxyz'),tf.transformations.euler_from_quaternion(q,'rzyx')
    print np.array(R.getQuaternion_matrix(q))
    print type(np.array(R.getQuaternion_matrix(q)))
    print(R.M)