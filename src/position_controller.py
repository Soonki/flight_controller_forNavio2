import numpy as np
from Robot_parameters import Drone_only as D
import tf
import matlab.engine


class position_controller():
    def __init__(self,D):
        self.D=D
        self.U=np.zeros((4,))
        self.Ts=0.01
        self.a0=(self.D.Ixx-self.D.Iyy)/self.D.Izz
        self.a1=(self.D.Izz-self.D.Ixx)/self.D.Iyy
        self.a2=(self.D.Iyy-self.D.Izz)/self.D.Ixx

        C10=600
        C11=25
        C12=25
        C20=15
        C21=3
        C22=3
        ram0=0.01
        ram1=0.01
        ram2=0.01

        self.K10=(1+ram0+C10*C20)
        self.K11=(1+ram1+C11*C21)
        self.K12=(1+ram2+C12*C22)
        self.K20=C10+C20
        self.K21=C11+C21
        self.K22=C12+C22
        self.K30=ram0*C20
        self.K31=ram1*C21
        self.K32=ram2*C22


        C13=2
        C14=2
        C15=10
        
        C23=0.5
        C24=0.5
        C25=1.5
        
        ram3=0.01
        ram4=0.01
        ram5=0.01

        self.K13=(1+ram3+C13*C23)
        self.K14=(1+ram4+C14*C24)
        self.K15=(1+ram5+C15*C25)
        self.K23=C13+C23
        self.K24=C14+C24
        self.K25=C15+C25
        self.K33=ram3*C23
        self.K34=ram4*C24
        self.K35=ram5*C25


        self.IMAX=1000
        self.IMIN=-1000

        self.ref=np.zeros((6,),dtype=float)#PHI,P
        self.p=np.zeros((3,))
        self.pos_int=np.zeros((3,))

        #Admittance filetr parameters
        self.M_E=np.array([[0.2,0.0,0.0],[0.0,0.2,0.0],[0.0,0.0,0.05]])
        self.D_E=np.array([[3.35,0.0,0.0],[0.0,3.35,0.0],[0.0,0.0,3.35]])
        self.K_E=np.array([[10.0,0.0,0.0],[0.0,10.0,0.0],[0.0,0.0,5.0]])

        #self.M_E=np.array([[1,0.0,0.0],[0.0,1,0.0],[0.0,0.0,1]])
        #self.D_E=np.array([[50,0.0,0.0],[0.0,50,0.0],[0.0,0.0,50]])
        #self.K_E=np.array([[80,0.0,0.0],[0.0,80,0.0],[0.0,0.0,80]])

        #self.M_E=np.array([[1,0.0,0.0],[0.0,1,0.0],[0.0,0.0,1]])
        #self.D_E=np.array([[2,0.0,0.0],[0.0,2,0.0],[0.0,0.0,2]])
        #self.K_E=np.array([[4,0.0,0.0],[0.0,4,0.0],[0.0,0.0,4]])

        #self.M_E=np.array([[4,0.0,0.0],[0.0,4,0.0],[0.0,0.0,4]])
        #self.D_E=np.array([[55,0.0,0.0],[0.0,55,0.0],[0.0,0.0,55]])
        #self.K_E=np.array([[80,0.0,0.0],[0.0,80,0.0],[0.0,0.0,80]])

        A=np.array([[self.M_E[0,0],self.M_E[0,1],self.M_E[0,2],self.D_E[0,0],self.D_E[0,1],self.D_E[0,2]],\
                   [self.M_E[1,0],self.M_E[1,1],self.M_E[1,2],self.D_E[1,0],self.D_E[1,1],self.D_E[1,2]],\
                   [self.M_E[2,0],self.M_E[2,1],self.M_E[2,2],self.D_E[2,0],self.D_E[2,1],self.D_E[2,2]],\
                   [0.0,0.0,0.0,1.0,0.0,0.0],[0.0,0.0,0.0,0.0,1.0,0.0],[0.0,0.0,0.0,0.0,0.0,1.0]])
        self.inv_A=np.linalg.inv(A)

        self.B=np.zeros((6,6))
        self.B[0:3,3:6]=self.K_E
        self.B[3:6,0:3]=-np.eye(3)

        self.v_dot=np.zeros((6,))
        self.v=np.zeros((6,))

    def setDrone(self,D):
        self.D=D

    def setReference(self,ref):
        self.ref=ref

    def updateIntegral(self):
        self.p=self.p+self.ref[0:3]-self.D.phi
        self.pos_int=self.pos_int+self.ref[3:6]-self.D.pos

        for i in range(0,2):
            if self.p[i]>self.IMAX:
                self.p[i]=self.IMAX
            elif self.p[i]<self.IMIN:
                self.p[i]=self.IMIN
            if self.pos_int[i]>self.IMAX:
                self.pos_int[i]=self.IMAX
            elif self.pos_int[i]<self.IMIN:
                self.pos_int[i]=self.IMIN

    def updateDesireTorque(self):
        self.U[1]=-self.D.Izz*(self.K10*(self.ref[0]-self.D.phi[0])-self.K20*self.D.phi_dot[0]+self.K30*self.p[0]-self.a0*self.D.phi_dot[1]*self.D.phi_dot[2])
        self.U[2]=-self.D.Iyy*(self.K11*(self.ref[1]-self.D.phi[1])-self.K21*self.D.phi_dot[1]+self.K31*self.p[1]-self.a1*self.D.phi_dot[2]*self.D.phi_dot[0])
        self.U[3]=self.D.Ixx*(self.K12*(self.ref[2]-self.D.phi[2])-self.K22*self.D.phi_dot[2]+self.K32*self.p[2]-self.a2*self.D.phi_dot[0]*self.D.phi_dot[1])

    def updateDesireTorque_link(self):
        self.U[1]=-self.D.Izz*(self.K10*(self.ref[0]-self.D.phi[0])-self.K20*self.D.phi_dot[0]+self.K30*self.p[0]-self.a0*self.D.phi_dot[1]*self.D.phi_dot[2])#+self.D.force[5]
        self.U[2]=-self.D.Iyy*(self.K11*(self.ref[1]-self.D.phi[1])-self.K21*self.D.phi_dot[1]+self.K31*self.p[1]-self.a1*self.D.phi_dot[2]*self.D.phi_dot[0])#+self.D.force[4]
        self.U[3]=self.D.Ixx*(self.K12*(self.ref[2]-self.D.phi[2])-self.K22*self.D.phi_dot[2]+self.K32*self.p[2]-self.a2*self.D.phi_dot[0]*self.D.phi_dot[1])\
                    -self.D.force[1]*(self.D.l_link0+(self.D.l_link1+self.D.l_ring)*np.cos(self.D.joint))-self.D.force[2]*(self.D.l_link1+self.D.l_ring)*np.sin(self.D.joint)\
                    -((self.D.m_link1+self.D.m_ring)*self.D.l_link1g+self.D.m_ring*self.D.l_ringg)*9.81*np.sin(self.D.joint)
                    #+self.D.force[3]

    def updateDesireThrustforAltitude(self):
        self.U[0]=self.D.M*(self.K15*(self.ref[5]-self.D.pos[2])+self.K25*(self.v[2]-self.D.pos_dot[2])+self.K35*self.pos_int[2]+9.81+self.v_dot[2])/self.D.Cth/self.D.Cps-self.D.force[2]

    def Mixing(self):
        self.Thrusts=np.dot(self.D.Mixing,self.U)
        return self.Thrusts

    def updateDesirePhiforPosition(self):
        self.Ux=self.D.M/self.U[0]*(self.K13*(self.ref[3]-self.D.pos[0])+self.K23*(self.v[0]-self.D.pos_dot[0])+self.K33*self.pos_int[0]+self.v_dot[0])+self.D.force[0]
        self.Uy=self.D.M/self.U[0]*(self.K14*(self.ref[4]-self.D.pos[1])+self.K24*(self.v[1]-self.D.pos_dot[1])+self.K34*self.pos_int[1]+self.v_dot[1])+self.D.force[1]
        if abs(self.Ux*np.sin(self.ref[0])-self.Uy*np.cos(self.ref[0]))>0.5:
            self.ref[2]=np.sign(self.Ux*np.sin(self.ref[0])-self.Uy*np.cos(self.ref[0]))*np.pi/6
        else:
            self.ref[2]=np.arcsin(self.Ux*np.sin(self.ref[0])-self.Uy*np.cos(self.ref[0]))

        if abs((self.Ux*np.cos(self.ref[0])+self.Uy*np.sin(self.ref[0]))/np.cos(self.ref[2]))>0.5:
            self.ref[1]=np.sign((self.Ux*np.cos(self.ref[0])+self.Uy*np.sin(self.ref[0]))/np.cos(self.ref[2]))*np.pi/6
        else:
            self.ref[1]=np.arcsin((self.Ux*np.cos(self.ref[0])+self.Uy*np.sin(self.ref[0]))/np.cos(self.ref[2]))

    def Attitudecontrol(self,D,ref):
        self.setReference(ref)
        self.setDrone(D)
        self.updateIntegral()
        self.updateDesireThrustforAltitude()
        self.updateDesireTorque()
        self.Mixing()

    def Positioncontrol(self,D,ref):
        self.setReference(ref)
        self.setDrone(D)
        self.updateIntegral()
        self.updateDesireThrustforAltitude()
        self.updateDesirePhiforPosition()
        self.updateDesireTorque()
        self.Mixing()

    def updateEndpos(self):
        self.R=tf.transformations.euler_matrix(self.D.phi[0],self.D.phi[1],self.D.phi[2],'rzyx')
        self.pos_end=self.D.pos+self.R[0:3,2]*self.D.l

    def updateEndpos_dot(self):
        pos_end_dot_r=np.array([-self.D.phi_dot[0]*self.R[1,2]+self.D.phi_dot[1]*self.D.Cph*self.D.Cth*self.D.Cps-self.D.phi_dot[2]*self.R[0,1],\
                                 self.D.phi_dot[0]*self.R[0,2]+self.D.phi_dot[1]*self.D.Sph*self.D.Cth*self.D.Cps-self.D.phi_dot[2]*self.R[1,1],\
                                 -self.D.phi_dot[1]*self.D.Sth*self.D.Cps-self.D.phi_dot[2]*self.D.Cth*self.D.Sps])
        self.pos_end_dot=self.D.pos_dot+pos_end_dot_r

    def overwriteforEnd(self):
        self.D.pos=self.pos_end
        #self.D.pos_dot=self.pos_end_dot

    def Positioncontrol_end(self,D,ref):
        self.setReference(ref)
        self.setDrone(D)
        self.updateEndpos()
        self.updateEndpos_dot()
        self.overwriteforEnd()
        self.updateIntegral()
        self.updateDesireThrustforAltitude()
        self.updateDesirePhiforPosition()
        self.updateDesireTorque()
        self.Mixing()

    def Admittance_filter(self,dt):
        z=np.zeros((6,))
        z[0:3]=np.dot(self.K_E,self.ref[3:6])-np.dot(self.D.R[0:3,0:3],self.D.force[0:3])
        self.v_dot=np.dot(self.inv_A,z-np.dot(self.B,self.v))
        self.v=self.v+self.v_dot*dt
        self.ref[3:6]=self.v[3:6]

    def Admittance_control_end(self,D,ref,dt):
        self.setReference(ref)
        self.Admittance_filter(dt)
        self.setDrone(D)
        self.updateEndpos()
        self.updateEndpos_dot()
        self.overwriteforEnd()
        self.updateIntegral()
        self.updateDesireThrustforAltitude()
        self.updateDesirePhiforPosition()
        self.updateDesireTorque()
        self.Mixing()

    def Admittance_control(self,D,ref,dt):
        self.setReference(ref)
        self.Admittance_filter(dt)
        self.setDrone(D)
        self.updateIntegral()
        self.updateDesireThrustforAltitude()
        self.updateDesirePhiforPosition()
        self.updateDesireTorque()
        self.Mixing()

    def Admittance_control_link(self,D,ref,dt):
        self.setReference(ref)
        self.Admittance_filter(dt)
        self.setDrone(D)
        self.updateIntegral()
        self.updateDesireThrustforAltitude()
        self.updateDesirePhiforPosition()
        self.updateDesireTorque_link()
        self.Mixing()

    def Position_control_link(self,D,ref):
        self.setReference(ref)
        self.setDrone(D)
        self.updateIntegral()
        self.updateDesireThrustforAltitude()
        self.updateDesirePhiforPosition()
        self.updateDesireTorque_link()
        self.Mixing()


class position_controller2():
    

    def __init__(self,D):
        self.D=D
        self.U=np.zeros((4,))
        self.Ts=0.01
        self.a0=(self.D.Ixx-self.D.Iyy)/self.D.Izz
        self.a1=(self.D.Izz-self.D.Ixx)/self.D.Iyy
        self.a2=(self.D.Iyy-self.D.Izz)/self.D.Ixx

        C10=600
        C11=45
        C12=45
        C20=15
        C21=3
        C22=3
        ram0=0.01
        ram1=0.01
        ram2=0.01

        self.K10=(1+ram0+C10*C20)
        self.K11=(1+ram1+C11*C21)
        self.K12=(1+ram2+C12*C22)
        self.K20=C10+C20
        self.K21=C11+C21
        self.K22=C12+C22
        self.K30=ram0*C20
        self.K31=ram1*C21
        self.K32=ram2*C22


        C13=4
        C14=4
        C15=10
        
        C23=0.5
        C24=0.5
        C25=1.5
        
        ram3=0.01
        ram4=0.01
        ram5=0.01

        self.K13=(1+ram3+C13*C23)
        self.K14=(1+ram4+C14*C24)
        self.K15=(1+ram5+C15*C25)
        self.K23=C13+C23
        self.K24=C14+C24
        self.K25=C15+C25
        self.K33=ram3*C23
        self.K34=ram4*C24
        self.K35=ram5*C25


        self.IMAX=1000
        self.IMIN=-1000

        self.ref=np.zeros((6,),dtype=float)#PHI,P
        self.p=np.zeros((3,))
        self.pos_int=np.zeros((3,))

        #Admittance filetr parameters
        self.M_E=np.array([[0.2,0.0,0.0],[0.0,0.2,0.0],[0.0,0.0,0.05]])
        self.D_E=np.array([[3.35,0.0,0.0],[0.0,3.35,0.0],[0.0,0.0,3.35]])
        self.K_E=np.array([[10.0,0.0,0.0],[0.0,10.0,0.0],[0.0,0.0,5.0]])

        #self.M_E=np.array([[1,0.0,0.0],[0.0,1,0.0],[0.0,0.0,1]])
        #self.D_E=np.array([[50,0.0,0.0],[0.0,50,0.0],[0.0,0.0,50]])
        #self.K_E=np.array([[80,0.0,0.0],[0.0,80,0.0],[0.0,0.0,80]])

        #self.M_E=np.array([[1,0.0,0.0],[0.0,1,0.0],[0.0,0.0,1]])
        #self.D_E=np.array([[2,0.0,0.0],[0.0,2,0.0],[0.0,0.0,2]])
        #self.K_E=np.array([[4,0.0,0.0],[0.0,4,0.0],[0.0,0.0,4]])

        #self.M_E=np.array([[4,0.0,0.0],[0.0,4,0.0],[0.0,0.0,4]])
        #self.D_E=np.array([[55,0.0,0.0],[0.0,55,0.0],[0.0,0.0,55]])
        #self.K_E=np.array([[80,0.0,0.0],[0.0,80,0.0],[0.0,0.0,80]])

        A=np.array([[self.M_E[0,0],self.M_E[0,1],self.M_E[0,2],self.D_E[0,0],self.D_E[0,1],self.D_E[0,2]],\
                   [self.M_E[1,0],self.M_E[1,1],self.M_E[1,2],self.D_E[1,0],self.D_E[1,1],self.D_E[1,2]],\
                   [self.M_E[2,0],self.M_E[2,1],self.M_E[2,2],self.D_E[2,0],self.D_E[2,1],self.D_E[2,2]],\
                   [0.0,0.0,0.0,1.0,0.0,0.0],[0.0,0.0,0.0,0.0,1.0,0.0],[0.0,0.0,0.0,0.0,0.0,1.0]])
        self.inv_A=np.linalg.inv(A)

        self.B=np.zeros((6,6))
        self.B[0:3,3:6]=self.K_E
        self.B[3:6,0:3]=-np.eye(3)

        self.v_dot=np.zeros((6,))
        self.v=np.zeros((6,))
        self.eng = matlab.engine.start_matlab()
        self.Thrusts=np.zeros((4,))

    def setDrone(self,D):
        self.D=D

    def getStateVector(self):
        return np.r_[self.D.pos[1],self.D.pos[0],-self.D.pos[2],-self.D.phi[0],self.D.phi[2],self.D.phi[1],\
                  self.D.pos_dot[1],self.D.pos_dot[0],-self.D.pos_dot[2],-self.D.phi_dot[0],self.D.phi_dot[2],self.D.phi_dot[1],\
                  np.pi,0.0,self.D.joint-np.pi/2,self.D.joint_dot]

    def getInputs(self):
        tau=self.Thrusts*0.01
        return np.r_[self.Thrusts,-tau[0],tau[1],-tau[2],tau[3],0,0]

    def setDynamics(self):
        x=self.getStateVector()
        self.Dynamics=self.eng.Dynamics(x.tolist())
        self.eng.workspace["D"]=self.Dynamics

    def updateDynamics(self):
        x=self.getStateVector()
        u=self.getInputs()
        self.eng.calculateMatrix(self.Dynamics,matlab.double(x.tolist()),matlab.double(u.tolist()),2,nargout=0)
        self.Mq=np.asarray(self.eng.eval("D.Mq"))
        self.Cq=np.asarray(self.eng.eval("D.Cq"))
        self.Gq=np.asarray(self.eng.eval("D.Gq"))

        #self.Mq[0,:],self.Mq[1,:],self.Mq[4,:],self.Mq[5,:]=self.Mq[1,:],self.Mq[0,:],self.Mq[5,:],self.Mq[4,:]
        #self.Mq[2,:],self.Mq[3,:]=-self.Mq[2,:],-self.Mq[3,:]

        #self.Cq[0],self.Cq[1],self.Cq[4],self.Cq[5]=self.Cq[1],self.Cq[0],self.Cq[5],self.Cq[4]
        #self.Cq[2],self.Cq[3]=-self.Cq[2],-self.Cq[3]

        #self.Gq[0],self.Gq[1],self.Gq[4],self.Gq[5]=self.Gq[1],self.Gq[0],self.Gq[5],self.Gq[4]
        #self.Gq[2],self.Gq[3]=-self.Gq[2],-self.Gq[3]

    def setReference(self,ref):
        self.ref=ref

    def updateIntegral(self):
        self.p=self.p+self.ref[0:3]-self.D.phi
        self.pos_int=self.pos_int+self.ref[3:6]-self.D.pos

        for i in range(0,2):
            if self.p[i]>self.IMAX:
                self.p[i]=self.IMAX
            elif self.p[i]<self.IMIN:
                self.p[i]=self.IMIN
            if self.pos_int[i]>self.IMAX:
                self.pos_int[i]=self.IMAX
            elif self.pos_int[i]<self.IMIN:
                self.pos_int[i]=self.IMIN

    def updateDesireTorque(self):
        self.U[1]=-self.D.Izz*(self.K10*(self.ref[0]-self.D.phi[0])-self.K20*self.D.phi_dot[0]+self.K30*self.p[0]-self.a0*self.D.phi_dot[1]*self.D.phi_dot[2])
        self.U[2]=-self.D.Iyy*(self.K11*(self.ref[1]-self.D.phi[1])-self.K21*self.D.phi_dot[1]+self.K31*self.p[1]-self.a1*self.D.phi_dot[2]*self.D.phi_dot[0])
        self.U[3]=self.D.Ixx*(self.K12*(self.ref[2]-self.D.phi[2])-self.K22*self.D.phi_dot[2]+self.K32*self.p[2]-self.a2*self.D.phi_dot[0]*self.D.phi_dot[1])

    def updateDesireTorque_link(self):
        self.U[1]=-self.Mq[3,3]*(self.K10*(self.ref[0]-self.D.phi[0])-self.K20*self.D.phi_dot[0]+self.K30*self.p[0])-self.Cq[3]-self.Gq[3]-self.D.tau[2]
        self.U[2]=-self.Mq[4,4]*(self.K11*(self.ref[1]-self.D.phi[1])-self.K21*self.D.phi_dot[1]+self.K31*self.p[1])-self.Cq[4]-self.Gq[4]-self.D.tau[1]
        self.U[3]=self.Mq[5,5]*(self.K12*(self.ref[2]-self.D.phi[2])-self.K22*self.D.phi_dot[2]+self.K32*self.p[2])+self.Cq[5]+self.Gq[5]+self.D.tau[0]
                    #-self.D.force[1]*(self.D.l_link0+(self.D.l_link1+self.D.l_ring)*np.cos(self.D.joint))-self.D.force[2]*(self.D.l_link1+self.D.l_ring)*np.sin(self.D.joint)\
                    #-((self.D.m_link1+self.D.m_ring)*self.D.l_link1g+self.D.m_ring*self.D.l_ringg)*9.81*np.sin(self.D.joint)
                    #+self.D.force[3]

    def updateDesireThrustforAltitude(self):
        self.U[0]=self.D.M*(self.K15*(self.ref[5]-self.D.pos[2])+self.K25*(self.v[2]-self.D.pos_dot[2])+self.K35*self.pos_int[2]+9.81+self.v_dot[2])/self.D.Cth/self.D.Cps+self.D.force[2]

    def Mixing(self):
        self.Thrusts=np.dot(self.D.Mixing,self.U)
        return self.Thrusts

    def updateDesirePhiforPosition(self):
        self.Ux=(self.Mq[0,0]*(self.K13*(self.ref[3]-self.D.pos[0])+self.K23*(self.v[0]-self.D.pos_dot[0])+self.K33*self.pos_int[0]+self.v_dot[0])+self.Cq[0]+self.Gq[0]+self.D.force[0])/self.U[0]
        self.Uy=(self.Mq[1,1]*(self.K14*(self.ref[4]-self.D.pos[1])+self.K24*(self.v[1]-self.D.pos_dot[1])+self.K34*self.pos_int[1]+self.v_dot[1])+self.Cq[1]+self.Gq[1]+self.D.force[1])/self.U[0]
        if abs(self.Ux*np.sin(self.ref[0])-self.Uy*np.cos(self.ref[0]))>0.5:
            self.ref[2]=np.sign(self.Ux*np.sin(self.ref[0])-self.Uy*np.cos(self.ref[0]))*np.pi/6
        else:
            self.ref[2]=np.arcsin(self.Ux*np.sin(self.ref[0])-self.Uy*np.cos(self.ref[0]))

        if abs((self.Ux*np.cos(self.ref[0])+self.Uy*np.sin(self.ref[0]))/np.cos(self.ref[2]))>0.5:
            self.ref[1]=np.sign((self.Ux*np.cos(self.ref[0])+self.Uy*np.sin(self.ref[0]))/np.cos(self.ref[2]))*np.pi/6
        else:
            self.ref[1]=np.arcsin((self.Ux*np.cos(self.ref[0])+self.Uy*np.sin(self.ref[0]))/np.cos(self.ref[2]))

    def Attitudecontrol(self,D,ref):
        self.setReference(ref)
        self.setDrone(D)
        self.updateIntegral()
        self.updateDesireThrustforAltitude()
        self.updateDesireTorque()
        self.Mixing()

    def Positioncontrol(self,D,ref):
        self.setReference(ref)
        self.setDrone(D)
        self.updateIntegral()
        self.updateDesireThrustforAltitude()
        self.updateDesirePhiforPosition()
        self.updateDesireTorque()
        self.Mixing()

    def updateEndpos(self):
        self.R=tf.transformations.euler_matrix(self.D.phi[0],self.D.phi[1],self.D.phi[2],'rzyx')
        self.pos_end=self.D.pos+self.R[0:3,2]*self.D.l

    def updateEndpos_dot(self):
        pos_end_dot_r=np.array([-self.D.phi_dot[0]*self.R[1,2]+self.D.phi_dot[1]*self.D.Cph*self.D.Cth*self.D.Cps-self.D.phi_dot[2]*self.R[0,1],\
                                 self.D.phi_dot[0]*self.R[0,2]+self.D.phi_dot[1]*self.D.Sph*self.D.Cth*self.D.Cps-self.D.phi_dot[2]*self.R[1,1],\
                                 -self.D.phi_dot[1]*self.D.Sth*self.D.Cps-self.D.phi_dot[2]*self.D.Cth*self.D.Sps])
        self.pos_end_dot=self.D.pos_dot+pos_end_dot_r

    def overwriteforEnd(self):
        self.D.pos=self.pos_end
        #self.D.pos_dot=self.pos_end_dot

    def Positioncontrol_end(self,D,ref):
        self.setReference(ref)
        self.setDrone(D)
        self.updateEndpos()
        self.updateEndpos_dot()
        self.overwriteforEnd()
        self.updateIntegral()
        self.updateDesireThrustforAltitude()
        self.updateDesirePhiforPosition()
        self.updateDesireTorque()
        self.Mixing()

    def Admittance_filter(self,dt):
        z=np.zeros((6,))
        z[0:3]=np.dot(self.K_E,self.ref[3:6])-self.D.force[0:3]
        self.v_dot=np.dot(self.inv_A,z-np.dot(self.B,self.v))
        self.v=self.v+self.v_dot*dt
        self.ref[3:6]=self.v[3:6]

    def Admittance_control_end(self,D,ref,dt):
        self.setReference(ref)
        self.Admittance_filter(dt)
        self.setDrone(D)
        self.updateEndpos()
        self.updateEndpos_dot()
        self.overwriteforEnd()
        self.updateIntegral()
        self.updateDesireThrustforAltitude()
        self.updateDesirePhiforPosition()
        self.updateDesireTorque()
        self.Mixing()

    def Admittance_control(self,D,ref,dt):
        self.setReference(ref)
        self.Admittance_filter(dt)
        self.setDrone(D)
        self.updateIntegral()
        self.updateDesireThrustforAltitude()
        self.updateDesirePhiforPosition()
        self.updateDesireTorque()
        self.Mixing()

    def Admittance_control_link(self,D,ref,dt):
        self.setReference(ref)
        self.Admittance_filter(dt)
        self.setDrone(D)
        self.updateIntegral()
        self.updateDynamics()
        self.updateDesireThrustforAltitude()
        self.updateDesirePhiforPosition()
        self.updateDesireTorque_link()
        self.Mixing()

    def Position_control_link(self,D,ref):
        self.setReference(ref)
        self.setDrone(D)
        self.updateIntegral()
        self.updateDynamics()
        self.updateDesireThrustforAltitude()
        self.updateDesirePhiforPosition()
        self.updateDesireTorque_link()
        self.Mixing()

if __name__ == '__main__':
    from Robot_parameters import Drone_one_joint
    import time
    ref=np.array([0.0,0.0,0.0,1.0,3.0,0.2])
    dt=0.01
    D=Drone_one_joint()
    D.force=np.ones((6,))
    a=position_controller2(D)
    a.setDynamics()
    i=0
    while(i<100):
        t=time.time()
        a.updateDynamics()
        print(time.time()-t)
        i=i+1