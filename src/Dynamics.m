classdef Dynamics < handle
    
    %-----------------初期位置---------------------------------------
    %            
    %            -------------------
    %                    |
    %       --->x        |
    %      |             |--------
    %      z              
    
 properties(Constant)
    %--------------------------------Body--------------------------------------------------------------------------
    m0=0.612; L=0.23; Ix=0.002496; Iy= 0.00252;
    Iz=0.004265;
    g=9.81;
    %--------------------------------Arm1--------------------------------------------------------------------------
    m1=0.007019; L1=0.05; l1=0.025; I1=diag([2.5*10^(-5),1*10^(-4),2.5*10^(-5)]);tau1_s=100000;
    %--------------------------------Arm2---------------------------------------------------------------------------
    m2=0.022401; L2=0.145; l2=0.079; I2=diag([2.5*10^(-5),1*10^(-4),2.5*10^(-5)]);tau2_s=100000;
 end
 
 properties(GetAccess = public, SetAccess = public)
    x;
    Mq;
    Cq;
    Gq;
    Fq;
 end
 
 methods
     %---------------------インスタンス---------------------------------
     function C=Dynamics(x)
         C.x=x;
         C.Mq=zeros(8);
         C.Cq=zeros(8,1);
         C.Gq=zeros(8,1);
         C.Fq=zeros(8,1);
     end
     %-------------------------行列計算---------------------------- 
     function calculateMatrix(C,x,u,LEVEL)
         %% x=[x y z ph th ps x_dot y_dot z_dot ph_dot th_dot ps_dot eta1 eta1_dot eta2 eta2_dot]
              %%[1 2  3 4   5  6    7      8      9      10       11     12      13     14      15      16    ]
         %% u=[T1 T2 T3 T4 tau1 tau2 tau3 tau4 tau_et1 tau_et2]

%-----------------------------------------準備するでぇ-----------------------------------------------
 m=20;
Cph=round(cos(x(4)),m);Sph=round(sin(x(4)),m);
Cth=round(cos(x(5)),m);Sth=round(sin(x(5)),m);
Cps=round(cos(x(6)),m);Sps=round(sin(x(6)),m);
C1=round(cos(x(13)),m);S1=round(sin(x(13)),m);
C12=round(cos(x(13)+x(15)-pi/2),m);S12=round(sin(x(13)+x(15)-pi/2),m); C12r=round(cos(x(13)+x(15)),m);S12r=round(sin(x(13)+x(15)),m);
ph_dot=x(10); th_dot=x(11); ps_dot=x(12); et1_dot=x(14);et2_dot=x(16);
m1l1=C.m1*C.l1; m1L1=C.m1*C.L1; R=Rotation(x(4),x(5),x(6)); R1=[C1 0 S1; 0 1 0; -S1 0 C1]; R2=[C12r 0 S12r; 0 1 0; -S12r 0 C12r];
m2l2=C.m2*C.l2; m2L1=C.m2*C.L1; p=[x(1) x(2) x(3)]'; pdot=[x(7) x(8) x(9)]';
    %-------------------------------------回転行列関係の微分----------------------------------------
            Rdot=[-ph_dot*Sph*Cth-th_dot*Cph*Sth -ph_dot*R(2,2)+th_dot*Cph*Cth*Sps+ps_dot*R(1,3) -ph_dot*R(2,3)+th_dot*Cph*Cth*Cps-ps_dot*R(1,2);
                  ph_dot*Cph*Cth-th_dot*Sph*Sth  ph_dot*R(1,2)+th_dot*Sph*Cth*Sps+ps_dot*R(2,3)  ph_dot*R(1,3)+th_dot*Sph*Cth*Cps-ps_dot*R(2,2);
                  -th_dot*Cth                    0                                               -th_dot*Sth*Cps-ps_dot*Cth*Sps                  ];
            
            Rdotdot=[-ph_dot^2*Cph*Cth+ph_dot*th_dot*(Cph*Cth+Sph*Sth)-th_dot^2*Cph*Cth 0 -ph_dot*Rdot(2,3)+th_dot*(-ph_dot*Sph*Cth*Cps-th_dot*Cph*Sth*Cps-ps_dot*Cph*Cth*Sps)-ps_dot*Rdot(1,2);
                     -ph_dot^2*Sph*Cth-2*ph_dot*th_dot*Cph*Sth-th_dot^2*Sph*Cth         0 ph_dot*Rdot(1,3)+th_dot*(ph_dot*Cph*Cth*Cps-th_dot*Sph*Sth*Cps-ps_dot*Sph*Cth*Sps)-ps_dot*Rdot(2,2);
                     th_dot^2*Sth                                                       0 -th_dot^2*Cth*Cps+2*th_dot*ps_dot*Sth*Sps-ps_dot^2*Cth*Cps];
            
            R_ph=[-Sph*Cth -Sph*Sth*Sps-Cph*Cps -Sph*Sth*Cps+Cph*Sps;
                  Cph*Cth  Cph*Sth*Sps-Sph*Cps  Cph*Sth*Cps+Sph*Sps;
                  0        0                    0                   ];
              
            R_th=[-Cph*Sth Cph*Cth*Sps Cph*Cth*Cps;
                  -Sph*Sth Sph*Cth*Sps Sph*Cth*Cps;
                  -Cth     -Sth*Sps    -Sth*Cps   ];
            
            R_ps=[0 Cph*Sth*Cps+Sph*Sps -Cph*Sth*Sps+Sph*Cps;
                  0 Sph*Sth*Cps-Cph*Sps -Sph*Sth*Sps-Cph*Cps;
                  0 Cth*Cps             -Cth*Sps            ];                
                 
            Rdot_ph=[-ph_dot*Cph*Cth+th_dot*Sph*Sth 0 -ph_dot*R_ph(2,3)-th_dot*Sph*Cth*Cps-ps_dot*R_ph(1,2);
                     -ph_dot*Sph*Cth-th_dot*Cph*Sth 0  ph_dot*R_ph(1,3)+th_dot*Cph*Cth*Cps-ps_dot*R_ph(2,2);
                     0                              0  0                                                  ];
            Rdot_th=[ph_dot*Sph*Sth-th_dot*Cph*Cth   0 -ph_dot*R_th(2,3)-th_dot*Cph*Sth*Cps-ps_dot*R_th(1,2);
                     -ph_dot*Cph*Sth-th_dot*Sph*Cth  0  ph_dot*R_th(1,3)-th_dot*Sph*Sth*Cps-ps_dot*R_th(2,2);
                     th_dot*Sth                      0 -th_dot*Cth*Cps+ps_dot*Sth*Sps                  ]; 
                 
            Rdot_ps=[0 0 -ph_dot*R_ps(2,3)-th_dot*Cph*Cth*Sps-ps_dot*R_ps(1,2);
                  0 0  ph_dot*R_ps(1,3)-th_dot*Sph*Cth*Sps-ps_dot*R_ps(2,2);
                  0 0 th_dot*Sth*Sps-ps_dot*Cth*Cps                  ];
              
            Rdot_ph_dot=[Sph*Cth -R(2,2) -R(2,3);
                         Cph*Cth R(1,2)  R(1,3);
                         0 0 0];
            Rdot_th_dot=[-Cph*Sth Cph*Cth*Sps Cph*Cth*Cps;
                         -Sph*Sth Sph*Cth*Sps Sph*Cth*Cps;
                         -Cth     0           -Sth*Cps];
            Rdot_ps_dot=[0 R(1,3) -R(1,2);
                         0 R(2,3) -R(2,2);
                         0 0      -Cth*Sps];
                     
            Rdotdot_ph_dot=[Rdot(2,1) -Rdot(2,2) -Rdot(2,3);
                            Rdot(1,1)  Rdot(1,2)  Rdot(1,3);
                            0          0          0        ];
                        
            Rdotdot_th_dot=[ph_dot*Sph*Sth-th_dot*Cph*Cth  0 -ph_dot*Sph*Cth*Cps-th_dot*Cph*Sth*Cps-ps_dot*Cph*Cth*Sps;
                            -ph_dot*Cph*Sth-th_dot*Sph*Cth 0 ph_dot*Cph*Cth*Cps-th_dot*Sph*Sth*Cps-ps_dot*Sph*Cth*Sps;
                            th_dot*Sth                     0 -th_dot*Cth*Cps+ps_dot*Sth*Sps                           ];
                        
            Rdotdot_ps_dot=[0  Rdot(1,3) -Rdot(1,2);
                            0  Rdot(2,3) -Rdot(2,2);
                            0          0 th_dot*Sth*Sps-ps_dot*Cth*Cps];

%------------------------------------------Body-------------------------------------------------------
%-----------------------------------------------------------------------------------------------------
        %% Inertia Matrix
        Mb=zeros(8,8);
        Mb(1,1)=C.m0;
        Mb(2,2)=C.m0;
        Mb(3,3)=C.m0;
        Mb(4,4)=Sth^2*C.Ix+Cth^2*(Sps^2*C.Iy+Cps^2*C.Iz);
        Mb(4,5)=Sps*Cps*Cth*(C.Iy-C.Iz);
        Mb(4,6)=-Sth*C.Ix;
        Mb(5,4)=Mb(4,5);
        Mb(5,5)=Cps^2*C.Iy+Sps^2*C.Iz;
        Mb(6,4)=Mb(4,6);
        Mb(6,6)=C.Ix;
        %% h Matrix
        Cb=zeros(8,1);
        Cb(4)=2*Sth*Cth*(C.Ix-Sps^2*C.Iy-Cps^2*C.Iz)*ph_dot*th_dot-Cth*(C.Ix-(Cps^2-Sps^2)*C.Iy+(Cps^2-Sps^2)*C.Iz)*th_dot*ps_dot+2*Cth^2*Cps*Sps*(C.Iy-C.Iz)*ps_dot*ph_dot-Sth*Sps*Cps*(C.Iy-C.Iz)*th_dot^2;
        Cb(5)=2*Sth*Cps*Sps*(-C.Iy+C.Iz)*ph_dot*th_dot+Cth*(C.Ix+Cps*(Cps-Sps)*C.Iy+(Sps^2-Cps^2)*C.Iz)*ph_dot*ps_dot+(Cps*(Cps-Sps)*C.Iy+2*Cps*Sps*C.Iz)*th_dot*ps_dot-Sth*(-C.Ix+Cth*(Sps^2*C.Iy-Cps^2*C.Iz))*ph_dot^2;
        Cb(6)=Cth*(-C.Ix-(Cps^2-Sps^2)*C.Iy+(Cps^2-Sps^2)*C.Iz)*ph_dot*th_dot-Cth^2*Sps*Cps*(C.Iy-C.Iz)*ph_dot^2+Cps*Sps*(C.Iy-C.Iz)*th_dot^2;
        %% Gravity Matrix
        Gb=zeros(8,1);
        Gb(3)=-C.m0*C.g;
        %% Force function
        %後で消す
        F=(u(1)+u(2)+u(3)+u(4));
        Fb(1:3,1)=R*[0;0;-F];
        Fb(4,1)=-u(5)+u(6)-u(7)+u(8);
        Fb(5,1)=(u(1)-u(2)-u(3)+u(4))*C.L/4*sqrt(2);
        Fb(6,1)=(-u(1)-u(2)+u(3)+u(4))*C.L/4*sqrt(2);
        Fb(7,1)=u(9);
        Fb(8,1)=u(10);
%------------------------------------------------------------------------------------------------------

%----------------------------------Arm-----------------------------------------------------------------
        %--------------------------1LINK-------------------------------------------------------------
            %----------------------並進エネルギについて-------------------------------------------
            p1b=[S1 0 C1]'.*C.l1;
            p1bdot=[C1 0 -S1]'.*et1_dot*C.l1;
            p1bdotdot=[-S1 0 -C1]'.*et1_dot^2*C.l1;
            p1bdot_et1_dot=[C1 0 -S1]'.*C.l1;
            p1bdotdot_et1_dot=-[S1 0 C1]'.*et1_dot*C.l1;
            p1bdot_et1=-[S1 0 C1]'.*et1_dot*C.l1;
            p1b_et1=[C1 0 -S1]'.*C.l1;
            
            p1dot=pdot+Rdot*p1b+R*p1bdot;
            p1dotdot=Rdotdot*p1b+2*Rdot*p1bdot+R*p1bdotdot;
            
            p1dotdot_x_dot=zeros(3,1);
            p1dotdot_y_dot=zeros(3,1);
            p1dotdot_z_dot=zeros(3,1);
            p1dotdot_ph_dot=Rdotdot_ph_dot*p1b+Rdot_ph_dot*p1bdot;
            p1dotdot_th_dot=Rdotdot_th_dot*p1b+Rdot_th_dot*p1bdot;
            p1dotdot_ps_dot=Rdotdot_ps_dot*p1b+Rdot_ps_dot*p1bdot;
            p1dotdot_et1_dot=R*p1bdotdot_et1_dot+Rdot*p1bdot_et1_dot;
            
            p1dot_x_dot=[1 0 0]';
            p1dot_y_dot=[0 1 0]';
            p1dot_z_dot=[0 0 1]';
            p1dot_ph_dot=Rdot_ph_dot*p1b;
            p1dot_th_dot=Rdot_th_dot*p1b;
            p1dot_ps_dot=Rdot_ps_dot*p1b;
            p1dot_et1_dot=R*p1bdot_et1_dot;
            
            p1dot_x=zeros(3,1);
            p1dot_y=zeros(3,1);
            p1dot_z=zeros(3,1);
            p1dot_ph=Rdot_ph*p1b+R_ph*p1bdot;
            p1dot_th=Rdot_th*p1b+R_th*p1bdot;
            p1dot_ps=Rdot_ps*p1b+R_ps*p1bdot;
            p1dot_et1=Rdot*p1b_et1+R*p1bdot_et1;
            
            %-----------------------慣性行列-------------------------------------------------------
            M1t=zeros(8);
            M1t(1,:)=[C.m1 0    0    (-Sph*Cth*S1-R(2,3)*C1)*m1l1 (-Cph*Sth*S1+Cph*Cth*Cps*C1)*m1l1 -R(1,2)*C1*m1l1  (Cph*Cth*C1-R(1,3)*S1)*m1l1   0];
            M1t(2,:)=[0    C.m1 0    (Cph*Cth*S1+R(1,3)*C1)*m1l1  (-Sph*Sth*S1+Sph*Cth*Cps*C1)*m1l1 -R(2,2)*C1*m1l1  (R(2,1)*C1-R(2,3)*S1)*m1l1 0];
            M1t(3,:)=[0    0    C.m1 0                            (-Cth*S1-Sth*Cps*C1)*m1l1         -Cth*Sps*C1*m1l1 (R(3,1)*C1-R(3,3)*S1)*m1l1 0];         
            M1t(4,:)=(Rdot_ph_dot*p1b)'*M1t(1:3,:);
            M1t(5,:)=(Rdot_th_dot*p1b)'*M1t(1:3,:);
            M1t(6,:)=(Rdot_ps_dot*p1b)'*M1t(1:3,:);
            M1t(7,:)=(R*p1bdot_et1_dot)'*M1t(1:3,:);
            %-----------------------コリオリ，遠心力--------------------------------------------------
            C1t=zeros(8,1);
            C1t(1,1)=C.m1*(p1dotdot_x_dot'*p1dot+p1dot_x_dot'*p1dotdot-p1dot_x'*p1dot);
            C1t(2,1)=C.m1*(p1dotdot_y_dot'*p1dot+p1dot_y_dot'*p1dotdot-p1dot_y'*p1dot);
            C1t(3,1)=C.m1*(p1dotdot_z_dot'*p1dot+p1dot_z_dot'*p1dotdot-p1dot_z'*p1dot);
            C1t(4,1)=C.m1*(p1dotdot_ph_dot'*p1dot+p1dot_ph_dot'*p1dotdot-p1dot_ph'*p1dot);
            C1t(5,1)=C.m1*(p1dotdot_th_dot'*p1dot+p1dot_th_dot'*p1dotdot-p1dot_th'*p1dot);
            C1t(6,1)=C.m1*(p1dotdot_ps_dot'*p1dot+p1dot_ps_dot'*p1dotdot-p1dot_ps'*p1dot);
            C1t(7,1)=C.m1*(p1dotdot_et1_dot'*p1dot+p1dot_et1_dot'*p1dotdot-p1dot_et1'*p1dot);
            
            %-------------------------回転エネルギについて-------------------------------------
            I1r=R1*C.I1*R1';
            I1rdot=[2*C1*S1*et1_dot*(-C.I1(1,1)+C.I1(3,3))    0 (C1^2-S1^2)*et1_dot*(-C.I1(1,1)+C.I1(3,3));
                    0                                         0 0                                        ;
                    (C1^2-S1^2)*et1_dot*(-C.I1(1,1)+C.I1(3,3)) 0 2*C1*S1*et1_dot*(C.I1(1,1)-C.I1(3,3))   ];
                
            I1r_et1=[2*C1*S1*(-C.I1(1,1)+C.I1(3,3))    0 (C1^2-S1^2)*(-C.I1(1,1)+C.I1(3,3));
                    0                                  0 0                                ;
                    (C1^2-S1^2)*(-C.I1(1,1)+C.I1(3,3))  0 2*C1*S1*(C.I1(1,1)-C.I1(3,3))   ];  
            
            om1=[-Sth*ph_dot+ps_dot;Cth*Sps*ph_dot+Cps*th_dot+et1_dot;Cth*Cps*ph_dot-Sps*th_dot];
            om1dot=[-th_dot*ph_dot*Cth; -ph_dot*th_dot*Sth*Sps-th_dot*ps_dot*Sps+ps_dot*ph_dot*Cth*Cps; -ph_dot*th_dot*Sth*Cps-th_dot*ps_dot*Cps-ps_dot*ph_dot*Cth*Sps];
            
            om1_ph_dot=[-Sth;Cth*Sps;Cth*Cps];
            om1_th_dot=[0;Cps;-Sps];
            om1_ps_dot=[1;0;0];
            om1_et1_dot=[0;1;0];
            
            om1dot_ph_dot=[-th_dot*Cth;-th_dot*Sth*Sps+ps_dot*Cth*Cps;-th_dot*Sth*Cps-ps_dot*Cth*Sps];
            om1dot_th_dot=[0;-ps_dot*Sps;-ps_dot*Cps];
            om1dot_ps_dot=zeros(3,1);
            om1dot_et1_dot=zeros(3,1);
            
            om1_ph=zeros(3,1);
            om1_th=[-Cth*ph_dot;-Sth*Sps*ph_dot;-Sth*Cps*ph_dot];
            om1_ps=[0;Cth*Cps*ph_dot-Sps*th_dot;-Cth*Sps*ph_dot-Cps*th_dot];
            om1_et1=zeros(3,1);
            
            
            %------------------------慣性行列------------------------------------------------
            M1r=zeros(8);
            M1r(4,:)=[0 0 0 I1r(1,1)*Sth^2-I1r(1,3)*Sth*Cth*Cps+I1r(2,2)*Cth^2*Sps^2-I1r(3,1)*Cth*Cps*Sth+I1r(3,3)*Cth^2*Cps^2 I1r(1,3)*Sth*Sps+I1r(2,2)*Cth*Sps*Cps-I1r(3,3)*Cth*Cps*Sps -I1r(1,1)*Sth+I1r(3,1)*Cth*Cps Cth*Sps*I1r(2,2) 0];
            M1r(5,:)=[0 0 0 I1r(2,2)*Cth*Sps*Cps-Sps*(-I1r(3,1)*Sth+I1r(3,3)*Cth*Cps) I1r(2,2)*Cps^2+I1r(3,3)*Sps^2 -I1r(3,1)*Sps I1r(2,2)*Cps 0];
            M1r(6,:)=[0 0 0 -Sth*I1r(1,1)+Cth*Cps*I1r(1,3) -I1r(1,3)*Sps I1r(1,1) 0 0];
            M1r(7,:)=[0 0 0 I1r(2,2)*Cth*Sps I1r(2,2)*Cps 0 I1r(2,2) 0];
            %-----------------------コリオリ，遠心力--------------------------------------------------
            C1r=zeros(8,1);
            C1r(4,1)=om1dot_ph_dot'*I1r*om1+om1_ph_dot'*I1rdot*om1+om1_ph_dot'*I1r*om1dot...
                     -0.5*(om1_ph'*I1r*om1+om1'*I1r*om1_ph);
            C1r(5,1)=om1dot_th_dot'*I1r*om1+om1_th_dot'*I1rdot*om1+om1_th_dot'*I1r*om1dot...
                     -0.5*(om1_th'*I1r*om1+om1'*I1r*om1_th);
            C1r(6,1)=om1dot_ps_dot'*I1r*om1+om1_ps_dot'*I1rdot*om1+om1_ps_dot'*I1r*om1dot...
                     -0.5*(om1_ps'*I1r*om1+om1'*I1r*om1_ps);
            C1r(7,1)=om1dot_et1_dot'*I1r*om1+om1_et1_dot'*I1rdot*om1+om1_et1_dot'*I1r*om1dot...
                     -0.5*(om1_et1'*I1r*om1+om1'*I1r_et1*om1+om1'*I1r*om1_et1);
            %--------------------------重力-------------------------------------------
            e3=[0 0 1]';
            G1=zeros(8,1);
            G1(3,1)=-C.m1*C.g;
            G1(4,1)=-C.m1*C.g*e3'*R_ph*p1b;
            G1(5,1)=-C.m1*C.g*e3'*R_th*p1b;
            G1(6,1)=-C.m1*C.g*e3'*R_ps*p1b;
            G1(7,1)=-C.m1*C.g*e3'*R*p1b_et1;
     %-------------------------------2LINK-----------------------------------------------
           %---------------------並進エネルギ--------------------------------------------
           p2b=[C.L1*S1-C.l2*S12;0;C.L1*C1-C.l2*C12];
           p2b_et1=[C.L1*C1-C.l2*C12;0;-C.L1*S1+C.l2*S12];
           p2b_et2=[-C.l2*C12;0;C.l2*S12];
           
           p2bdot=[C.L1*et1_dot*C1-C.l2*(et1_dot+et2_dot)*C12;0;-C.L1*et1_dot*S1+C.l2*(et1_dot+et2_dot)*S12];
           p2bdotdot=[-C.L1*et1_dot^2*S1+C.l2*(et1_dot+et2_dot)^2*S12;0;-C.L1*et1_dot^2*C1+C.l2*(et1_dot+et2_dot)^2*C12];
           
           p2bdot_et1=[-C.L1*et1_dot*S1+C.l2*(et1_dot+et2_dot)*S12;0;-C.L1*et1_dot*C1+C.l2*(et1_dot+et2_dot)*C12];
           p2bdot_et2=[C.l2*(et1_dot+et2_dot)*S12;0;C.l2*(et1_dot+et2_dot)*C12];
           
           p2bdot_et1_dot=[C.L1*C1-C.l2*C12;0;-C.L1*S1+C.l2*S12];
           p2bdot_et2_dot=[-C.l2*C12;0;C.l2*S12];
           p2bdotdot_et1_dot=[-C.L1*S1*et1_dot+C.l2*S12*(et1_dot+et2_dot);0;-C.L1*C1*et1_dot+C.l2*C12*(et1_dot+et2_dot)];
           p2bdotdot_et2_dot=[C.l2*S12*(et1_dot+et2_dot);0;C.l2*C12*(et1_dot+et2_dot)];
           
           p2dot=pdot+Rdot*p2b+R*p2bdot;
           p2dotdot=Rdotdot*p2b+2*Rdot*p2bdot+R*p2bdotdot;
           
           p2dot_x=zeros(3,1);
           p2dot_y=zeros(3,1);
           p2dot_z=zeros(3,1);
           p2dot_ph=Rdot_ph*p2b+R_ph*p2bdot;
           p2dot_th=Rdot_th*p2b+R_th*p2bdot;
           p2dot_ps=Rdot_ps*p2b+R_ps*p2bdot;
           p2dot_et1=R*p2bdot_et1+Rdot*p2b_et1;
           p2dot_et2=R*p2bdot_et2+Rdot*p2b_et2;
           
           p2dot_x_dot=[1;0;0];
           p2dot_y_dot=[0;1;0];
           p2dot_z_dot=[0;0;1];
           p2dot_ph_dot=Rdot_ph_dot*p2b;
           p2dot_th_dot=Rdot_th_dot*p2b;
           p2dot_ps_dot=Rdot_ps_dot*p2b;
           p2dot_et1_dot=R*p2bdot_et1_dot;
           p2dot_et2_dot=R*p2bdot_et2_dot;
           
           p2dotdot_x_dot=zeros(3,1);
           p2dotdot_y_dot=zeros(3,1);
           p2dotdot_z_dot=zeros(3,1);
           p2dotdot_ph_dot=Rdotdot_ph_dot*p2b+Rdot_ph_dot*p2bdot;
           p2dotdot_th_dot=Rdotdot_th_dot*p2b+Rdot_th_dot*p2bdot;
           p2dotdot_ps_dot=Rdotdot_ps_dot*p2b+Rdot_ps_dot*p2bdot;
           p2dotdot_et1_dot=Rdot*p2bdot_et1_dot+R*p2bdotdot_et1_dot;
           p2dotdot_et2_dot=Rdot*p2bdot_et2_dot+R*p2bdotdot_et2_dot;

           %-------------------慣性行列-------------------------------------------------
           M2t=zeros(8);
           M2t(1,:)=[1 0 0 -(C.L1*S1-C.l2*S12)*Sph*Cth-(C.L1*C1-C.l2*C12)*R(2,3) -(C.L1*S1-C.l2*S12)*Cph*Sth+(C.L1*C1-C.l2*C12)*Cph*Cth*Cps -(C.L1*C1-C.l2*C12)*R(1,2)  Cph*Cth*(C.L1*C1-C.l2*C12)+R(1,3)*(-C.L1*S1+C.l2*S12)   (-Cph*Cth*C12+R(1,3)*S12)*C.l2  ].*C.m2;
           M2t(2,:)=[0 1 0 (C.L1*S1-C.l2*S12)*Cph*Cth+(C.L1*C1-C.l2*C12)*R(1,3)  -(C.L1*S1-C.l2*S12)*Sph*Sth+(C.L1*C1-C.l2*C12)*Sph*Cth*Cps -(C.L1*C1-C.l2*C12)*R(2,2)  R(2,1)*(C.L1*C1-C.l2*C12)+R(2,3)*(-C.L1*S1+C.l2*S12) (-R(2,1)*C12+R(2,3)*S12)*C.l2].*C.m2;
           M2t(3,:)=[0 0 1 0                                                     -Cth*(C.L1*S1-C.l2*S12)-Sth*Cps*(C.L1*C1-C.l2*C12)         -Cth*Sps*(C.L1*C1-C.l2*C12) R(3,1)*(C.L1*C1-C.l2*C12)+R(3,3)*(-C.L1*S1+C.l2*S12) (-R(3,1)*C12+R(3,3)*S12)*C.l2].*C.m2;         
           M2t(4,:)=(Rdot_ph_dot*p2b)'*M2t(1:3,:);
           M2t(5,:)=(Rdot_th_dot*p2b)'*M2t(1:3,:);
           M2t(6,:)=(Rdot_ps_dot*p2b)'*M2t(1:3,:);
           M2t(7,:)=(R*p2bdot_et1_dot)'*M2t(1:3,:);
           M2t(8,:)=(R*p2bdot_et2_dot)'*M2t(1:3,:);
           %-----------------------コリオリ，遠心力---------------------------------------
           C2t=zeros(8,1);
           C2t(1,1)=C.m2*(p2dotdot_x_dot'*p2dot+p2dot_x_dot'*p2dotdot-p2dot_x'*p2dot);
           C2t(2,1)=C.m2*(p2dotdot_y_dot'*p2dot+p2dot_y_dot'*p2dotdot-p2dot_y'*p2dot);
           C2t(3,1)=C.m2*(p2dotdot_z_dot'*p2dot+p2dot_z_dot'*p2dotdot-p2dot_z'*p2dot);
           C2t(4,1)=C.m2*(p2dotdot_ph_dot'*p2dot+p2dot_ph_dot'*p2dotdot-p2dot_ph'*p2dot);
           C2t(5,1)=C.m2*(p2dotdot_th_dot'*p2dot+p2dot_th_dot'*p2dotdot-p2dot_th'*p2dot);
           C2t(6,1)=C.m2*(p2dotdot_ps_dot'*p2dot+p2dot_ps_dot'*p2dotdot-p2dot_ps'*p2dot);
           C2t(7,1)=C.m2*(p2dotdot_et1_dot'*p2dot+p2dot_et1_dot'*p2dotdot-p2dot_et1'*p2dot);
           C2t(8,1)=C.m2*(p2dotdot_et2_dot'*p2dot+p2dot_et2_dot'*p2dotdot-p2dot_et2'*p2dot);

           %----------------------回転エネルギ----------------------------------------------
            I2r=[C12r^2*C.I2(1,1)+S12r^2*C.I2(3,3) 0         -S12r*C12r*(C.I2(1,1)-C.I2(3,3));
                 0                                 C.I2(2,2) 0;
                 -S12r*C12r*(C.I2(1,1)-C.I2(3,3))  0         C12r^2*C.I2(1,1)+S12r^2*C.I2(3,3)];
            %I2r=R2*C.I2*R2';
                
            I2rdot=[2*C12r*S12r*(et1_dot+et2_dot)*(-C.I2(1,1)+C.I2(3,3))    0 (C12r^2-S12r^2)*(et1_dot+et2_dot)*(-C.I2(1,1)+C.I2(3,3));
                    0                                                       0               0                                        ;
                    (C12r^2-S12r^2)*(et1_dot+et2_dot)*(-C.I2(1,1)+C.I2(3,3)) 0  2*C12r*S12r*(et1_dot+et2_dot)*(C.I2(1,1)-C.I2(3,3))   ];
            
            I2r_et1=[2*C12r*S12r*(-C.I2(1,1)+C.I2(3,3))    0 (C12r^2-S12r^2)*(-C.I2(1,1)+C.I2(3,3));
                    0                                      0               0                      ;
                    (C12r^2-S12r^2)*(-C.I2(1,1)+C.I2(3,3))  0  2*C12r*S12r*(C.I2(1,1)-C.I2(3,3))  ];
            I2r_et2=I2r_et1; 
           
            om2=[-Sth*ph_dot+ps_dot;Cth*Sps*ph_dot+Cps*th_dot+et1_dot+et2_dot;Cth*Cps*ph_dot-Sps*th_dot];
            
            om2_ph_dot=[-Sth;Cth*Sps;Cth*Cps];
            om2_th_dot=[0;Cps;-Sps];
            om2_ps_dot=[1;0;0];
            om2_et1_dot=[0;1;0];
            om2_et2_dot=[0;1;0];
            
            om2dot_ph_dot=[-th_dot*Cth;-th_dot*Sth*Sps+ps_dot*Cth*Cps;-th_dot*Sth*Cps-ps_dot*Cth*Sps];
            om2dot_th_dot=[0;-ps_dot*Sps;-ps_dot*Cps];
            om2dot_ps_dot=zeros(3,1);
            om2dot_et1_dot=zeros(3,1);
            om2dot_et2_dot=zeros(3,1);
            
            om2_ph=zeros(3,1);
            om2_th=[-Cth*ph_dot;-Sth*Sps*ph_dot;-Sth*Cps*ph_dot];
            om2_ps=[0;Cth*Cps*ph_dot-Sps*th_dot;-Cth*Sps*ph_dot-Cps*th_dot];
            om2_et1=zeros(3,1);
            om2_et2=zeros(3,1);
            
            om2dot=[-th_dot*ph_dot*Cth; -ph_dot*th_dot*Sth*Sps-th_dot*ps_dot*Sps+ps_dot*ph_dot*Cth*Cps; -ph_dot*th_dot*Sth*Cps-th_dot*ps_dot*Cps-ps_dot*ph_dot*Cth*Sps];
           %------------------------慣性行列------------------------------------------------
            M2r=zeros(8);
            M2r(4,:)=[0 0 0 I2r(1,1)*Sth^2-I2r(1,3)*Sth*Cth*Cps+I2r(2,2)*Cth^2*Sps^2-I2r(3,1)*Sth*Cth*Cps+I2r(3,3)*Cth^2*Cps^2     I2r(1,3)*Sth*Sps+I2r(2,2)*Cth*Sps*Cps-I2r(3,3)*Cth*Cps*Sps -I2r(1,1)*Sth+I2r(3,1)*Cth*Cps Cth*Sps*I2r(2,2) Cth*Sps*I2r(2,2)];
            M2r(5,:)=[0 0 0 I2r(2,2)*Cth*Sps*Cps-Sps*(-I2r(3,1)*Sth+I2r(3,3)*Cth*Cps)                                              I2r(2,2)*Cps^2+I2r(3,3)*Sps^2                              -I2r(3,1)*Sps                  I2r(2,2)*Cps     I2r(2,2)*Cps    ];
            M2r(6,:)=[0 0 0 -Sth*I2r(1,1)+Cth*Cps*I2r(1,3)                                                                         -I2r(1,3)*Sps                                               I2r(1,1)                      0                0];
            M2r(7,:)=[0 0 0 I2r(2,2)*Cth*Sps                                                                                       I2r(2,2)*Cps                                                0                             I2r(2,2)         I2r(2,2)];
            M2r(8,:)=[0 0 0 I2r(2,2)*Cth*Sps                                                                                       I2r(2,2)*Cps                                                0                             I2r(2,2)         I2r(2,2)];
            %-----------------------コリオリ，遠心力--------------------------------------------------
            C2r=zeros(8,1);
            C2r(4,1)=om2dot_ph_dot'*I2r*om2+om2_ph_dot'*I2rdot*om2+om2_ph_dot'*I2r*om2dot...
                     -0.5*(om2_ph'*I2r*om2+om2'*I2r*om2_ph);
            C2r(5,1)=om2dot_th_dot'*I2r*om2+om2_th_dot'*I2rdot*om2+om2_th_dot'*I2r*om2dot...
                     -0.5*(om2_th'*I2r*om2+om2'*I2r*om2_th);
            C2r(6,1)=om2dot_ps_dot'*I2r*om2+om2_ps_dot'*I2rdot*om2+om2_ps_dot'*I2r*om2dot...
                     -0.5*(om2_ps'*I2r*om2+om2'*I2r*om2_ps);
            C2r(7,1)=om2dot_et1_dot'*I2r*om2+om2_et1_dot'*I2rdot*om2+om2_et1_dot'*I2r*om2dot...
                     -0.5*(om2_et1'*I2r*om2+om2'*I2r_et1*om2+om2'*I2r*om2_et1);
            C2r(8,1)=om2dot_et2_dot'*I2r*om2+om2_et2_dot'*I2rdot*om2+om2_et2_dot'*I2r*om2dot...
                     -0.5*(om2_et2'*I2r*om2+om2'*I2r_et2*om2+om2'*I2r*om2_et2);
           
            %--------------------------重力-------------------------------------------
            G2=zeros(8,1);
            G2(3,1)=-C.m2*C.g;
            G2(4,1)=-C.m2*C.g*e3'*R_ph*p2b;
            G2(5,1)=-C.m2*C.g*e3'*R_th*p2b;
            G2(6,1)=-C.m2*C.g*e3'*R_ps*p2b;
            G2(7,1)=-C.m2*C.g*e3'*R*p2b_et1;
            G2(8,1)=-C.m2*C.g*e3'*R*p2b_et2;
           %-----------------------合成----------------------------------------------------
            if LEVEL==0
                C.Mq=zeros(6,6);
                C.Cq=zeros(6,1);
                C.Gq=zeros(6,1);
                C.Fq=zeros(6,1);
                C.Mq=Mb(1:6,1:6);
                C.Cq=Cb(1:6,1);
                C.Gq=Gb(1:6,1);
                C.Fq=Fb(1:6,1);
            elseif LEVEL==1
                C.Mq=zeros(7);
                C.Cq=zeros(7,1);
                C.Gq=zeros(7,1);
                C.Fq=zeros(7,1);
                C.Mq=Mb(1:7,1:7)+M1t(1:7,1:7)+M1r(1:7,1:7);
                C.Cq=Cb(1:7,1)+C1t(1:7,1)+C1r(1:7,1);
                C.Gq=Gb(1:7,1)+G1(1:7,1);
                C.Fq=Fb(1:7,1);
                
                %-------------静止摩擦力----------------------
%                 F_all=C.Mq\(C.Fq-C.Cq-C.Gq);
%                 if F_all(7)<C.tau1_s && et1_dot==0
%                     C.Fq(7)=C.Fq(7)-(C.Fq(7)-C.Cq(7)-C.Gq(7));
%                 end
            elseif LEVEL==2
                C.Mq=zeros(8);
                C.Cq=zeros(8,1);
                C.Gq=zeros(8,1);
                C.Fq=zeros(8,1);
                C.Mq=Mb+M1t+M1r+M2t+M2r;
                C.Cq=Cb+C1t+C1r+C2t+C2r;
                C.Gq=Gb+G1+G2;
                C.Fq=Fb;
                %-------------静止摩擦力----------------------
%                 F_all=C.Mq\(C.Fq-C.Cq-C.Gq);
%                 if abs(F_all(7))<C.tau1_s && et1_dot==0
%                     C.Fq(7)=C.Fq(7)-(C.Fq(7)-C.Cq(7)-C.Gq(7));
%                     disp("ほい1");
%                     et1_dot
%                 end
%                 if abs(F_all(8))<C.tau1_s && et2_dot==0
%                     C.Fq(8)=C.Fq(8)-(C.Fq(8)-C.Cq(8)-C.Gq(8));
%                     disp("ほい2");
%                     et2_dot
%                 end
            end
     end
    %------------------------------------------------------------------
    %------------------------動く奴------------------------------------
    function Move(C,u,dt,LEVEL)
        %x=[x y z ph th ps x_dot y_dot z_dot ph_dot th_dot ps_dot eta1 eta1_dot et2 eta2_dot]
        calculateMatrix(C,C.x,u,LEVEL);
        if LEVEL==0
            C.x(13:16,:)=[];        
        elseif LEVEL==1
            C.x(15:16,:)=[];
            y=C.x;
            C.x(7)=y(13);
            C.x(8)=y(7);
            C.x(9)=y(8);
            C.x(10)=y(9);
            C.x(11)=y(10);
            C.x(12)=y(11);
            C.x(13)=y(12);         
        elseif LEVEL==2
            y=C.x;
            C.x(7)=y(13);
            C.x(8)=y(15);
            C.x(9)=y(7);
            C.x(10)=y(8);
            C.x(11)=y(9);
            C.x(12)=y(10);
            C.x(13)=y(11);
            C.x(14)=y(12);
            C.x(15)=y(14);
        end
        %x=[x y z ph th ps eta1 eta2 x_dot y_dot z_dot ph_dot th_dot ps_dot eta1_dot eta2_dot]
        num=length(C.x);
        A=zeros(num);
        A(1:num/2,1:num/2)=eye(num/2);
        A(num/2+1:num,num/2+1:num)=C.Mq;
        B=zeros(num,1);
        B(1:num/2)=-C.x(num/2+1:num);
        B(num/2+1:num)=C.Cq+C.Gq;
        Cc=zeros(num,1);
        Cc(num/2+1:num)=C.Fq;
        %disp(A)
        if isnan(det(A))
            disp("発散!!");
            disp(A)
            pause
        end
        C.x=C.x+A\(-B+Cc).*dt; 
        %C.x=round(C.x,7);
        %C.x(1:12,1)=zeros(12,1);
        %x=[x y z ph th ps x_dot y_dot z_dot ph_dot th_dot ps_dot eta1 eta1_dot eta1 eta2_dot]
        if LEVEL==0
            C.x(13:16,:)=zeros(4,1);
        elseif LEVEL==1
            y=C.x;
            C.x(13)=y(7);
            C.x(7)=y(8);
            C.x(8)=y(9);
            C.x(9)=y(10);
            C.x(10)=y(11);
            C.x(11)=y(12);
            C.x(12)=y(13);
            C.x(15:16)=zeros(2,1);
        elseif LEVEL==2
            y=C.x;
            C.x(15)=y(8);
            C.x(13)=y(7);
            C.x(7)=y(9);
            C.x(8)=y(10);
            C.x(9)=y(11);
            C.x(10)=y(12);
            C.x(11)=y(13);
            C.x(12)=y(14);
            C.x(14)=y(15);
        end
    end
 end
end