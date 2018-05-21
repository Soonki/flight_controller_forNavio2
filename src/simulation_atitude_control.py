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

Inputs=Float64MultiArray()
ref_r=Point()
log=Drone_log()
#-------------------------------------
#-----------Controller------------------
Ts=0.01
p=np.zeros((3,))
U=np.zeros((4,))


#-------------------------------------
#=============================================================================================
ref=np.array([0.0,0.0,0.0,0.0,0.0,0.845])
Cont.v[3:6]=ref[3:6]

Cont.setDynamics()
print('OK')
def callback(msg,p):
    ref[3:6]=np.array([msg.ref.x,msg.ref.y,msg.ref.z])
    Robot.update(msg)
    log.ref.x=ref[3]
    log.ref.y=ref[4]
    log.ref.z=ref[5]
    #Cont.Attitudecontrol(Robot,ref)
    #Cont.Positioncontrol(Robot,ref)
    #Cont.Positioncontrol_end(Robot,ref)
    #Cont.Admittance_control_end(Robot,ref,Ts)
    #Cont.Admittance_control(Robot,ref,Ts)
    Cont.Admittance_control_link(Robot,ref,Ts)
    #Cont.Position_control_link(Robot,ref)
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
    #print(Robot.phi)

if __name__ == '__main__':
    sub=rospy.Subscriber('info',Droneinfo_link,callback,p)
    rospy.spin()