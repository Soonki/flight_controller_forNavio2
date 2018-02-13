#!/usr/bin/env python
import numpy as np

class Drone_only():
    def __init__(self):
        self.M=0.612#[kg]
        self.L=0.29#[m]
        self.Lfh=0.09
        self.Lrh=0.10
        self.Lfv=0.115
        self.Lrv=0.107
        self.Ixx=0.0001#[kgm^2]
        self.Iyy=0.0001
        self.Izz=0.0002
        self.C=10
        self.Mixing=np.array([[      1,     1,      1,     1],\
                              [-self.C,self.C,-self.C,self.C],\
                              [self.Lfh,self.Lfh,-self.Lrh,-self.Lrh],\
                              [self.Lfv,-self.Lfv,-self.Lrv,self.Lrv]])
        self.Mixing=np.linalg.inv(self.Mixing)