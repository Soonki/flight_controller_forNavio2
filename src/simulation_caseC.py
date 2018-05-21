#!/usr/bin/env python
import time
import argparse 
import sys
import numpy as np
import rospy
import flight_functions as fcn
from std_msgs.msg import Float64MultiArray
from Robot_parameters import Drone_one_joint
from original_msgs.msg import Droneinfo_link,Drone_log
from geometry_msgs.msg import Point
import tf
from position_controller import position_controller2

Robot=Drone_one_joint()
Cont=position_controller2(Robot)

#=================================Initial setup==============================================
#-------------ROS-----------------------
FRQ=250.0
dt=1/FRQ

rospy.init_node('simulation_atitude_control')
pub1=rospy.Publisher('inputs',Float64MultiArray,queue_size=50)
pub2=rospy.Publisher('logging',Drone_log,queue_size=50)
pub3=rospy.Publisher('ref_r',Point,queue_size=50)
pub4=rospy.Publisher('ref',Point,queue_size=50)

Inputs=Float64MultiArray()
ref_r=Point()
ref_t=Point()
log=Drone_log()
#-------------------------------------
#-----------Controller------------------
Ts=0.01
p=np.zeros((3,))
U=np.zeros((4,))


#-------------------------------------
#=============================================================================================
ref=np.array([0.0,0.0,0.0,0.0,0.0,0.818])
Cont.v[3:6]=ref[3:6]

Cont.setDynamics()


class checker():

    def __init__(self):
        self.count=0

    def checkShock(self,R,ref_t,msg):
        if (self.count==0 and msg.joint.effort[0]>0.7):
            if (R.force[0] > 0.2):
                ref_t.x=0.65
                self.count=1
        elif (self.count==1 and msg.joint.effort[0]>0.7):
            ref_t.x=0.65
        return self.count

c=checker()
count=0
print('OK')
def callback(msg):
    ref_t.x=1.0
    ref_t.y=0.0
    ref_t.z=0.818
    ref[3:6]=np.array([msg.ref.x,msg.ref.y,msg.ref.z])
    Robot.update(msg)
    c.checkShock(Robot,ref_t,msg)
    log.ref.x=msg.ref.x
    log.ref.y=msg.ref.y
    log.ref.z=msg.ref.z
    #Cont.Attitudecontrol(Robot,ref)
    #Cont.Positioncontrol(Robot,ref)
    #Cont.Positioncontrol_end(Robot,ref)
    #Cont.Admittance_control_end(Robot,ref,Ts)
    #Cont.Admittance_control(Robot,ref,Ts)
    #Cont.Admittance_control_link(Robot,ref,Ts)
    Cont.Position_control_link(Robot,ref)
    Inputs.data=Cont.Thrusts
    log.odom.pose.pose.position.x=Robot.pos[0]
    log.odom.pose.pose.position.y=Robot.pos[1]
    log.odom.pose.pose.position.z=Robot.pos[2]
    log.odom.pose.pose.orientation.x=Robot.quaternion[0]
    log.odom.pose.pose.orientation.y=Robot.quaternion[1]
    log.odom.pose.pose.orientation.z=Robot.quaternion[2]
    log.odom.pose.pose.orientation.w=Robot.quaternion[3]
    log.odom.twist.twist.linear.x=Robot.pos_dot[0]
    log.odom.twist.twist.linear.y=Robot.pos_dot[1]
    log.odom.twist.twist.linear.z=Robot.pos_dot[2]
    log.odom.twist.twist.angular.x=Robot.phi_dot[2]
    log.odom.twist.twist.angular.y=Robot.phi_dot[1]
    log.odom.twist.twist.angular.z=Robot.phi_dot[0]
    log.load.force.x=Robot.force[0]
    log.load.force.y=Robot.force[1]
    log.load.force.z=Robot.force[2]
    log.joint.position=[Robot.joint]
    log.joint.velocity=[Robot.joint_dot]
    log.target_pos.x=Cont.ref[3]
    log.target_pos.y=Cont.ref[4]
    log.target_pos.z=Cont.ref[5]
    log.target_ang.x=Cont.ref[2]
    log.target_ang.y=Cont.ref[1]
    log.target_ang.z=Cont.ref[0]
    log.euler.x=Robot.phi[2]
    log.euler.y=Robot.phi[1]
    log.euler.z=Robot.phi[0]
    log.inputs.data=Cont.Thrusts
    ref_r.x=Cont.ref[3]
    ref_r.y=Cont.ref[4]
    ref_r.z=Cont.ref[5]
    pub1.publish(Inputs)
    pub2.publish(log)
    pub3.publish(ref_r)
    pub4.publish(ref_t)
    #print(Robot.phi)

if __name__ == '__main__':
    sub=rospy.Subscriber('info',Droneinfo_link,callback)
    rospy.spin()