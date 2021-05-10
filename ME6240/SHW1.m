% ME 6240 Short Homework #1
% Michael Swafford
% 02/05/2021
clear all
clc
%The goal for this short homework is to review some linear algebra conecpts
%in the context of a review of 3D rigid body dynamics

%Let’s say the orientation of some rigid body is defined as 3 subsequent 
%intrinsicEuler anglerotations, (𝜓,𝜃,𝜙), about the body’s (𝑧,𝑦,𝑥)axes, respectively,to 
%align it’s axeswith someinertial frame’s(𝑋,𝑌,𝑍)axes.The corresponding rotation matrix, 
%𝑹𝑊/𝐵,can be constructedby multiplying the intrinsicEuler angle rotation 
%matrices(below)in the appropriate order

%% 1 
%Provide RW/B as a function of (𝜓,𝜃,𝜙) by multiplying Euler angle roation matrices
%in the appropriate order
%use intrinsic rotation for rotating the frame w/respect to itself
%use extrinsic rotation for rotating the frame w/respect to an interial
%frame
disp("----------------------------- Problem 1 -----------------------------")
syms x y z
syms phi theta psi

rotx = [1 0 0;0 cos(phi) -sin(phi);0 sin(phi) cos(phi)];
roty = [cos(theta) 0 sin(theta);0 1 0;-sin(theta) 0 cos(theta)];
rotz = [cos(psi) -sin(psi) 0;sin(psi) cos(psi) 0;0 0 1];

B = [z y x]';

RWB = rotz*roty*rotx;
fprintf("The rotation matrix function to intrinsically rotate a rigid body (z,y,x)\n")
fprintf("by (psi,theta,phi) is:\n RWB = rotz*roty*rotx \n")
disp(RWB)
fprintf("\n")

%% 2 & 3 
disp("----------------------------- Problems 2 & 3 -----------------------------")

disp("2. S is a square matrix whoes transpose is equal to it's opposite")
disp("This must be the case since S + S' = 0")

RxT = [1 0 0;0 cos(phi) sin(phi);0 -sin(phi) cos(phi)];
drRx = diff(rotx,phi)
S = drRx*RxT

disp("3. For any angle phi S will evaluate to:")
s=[0 0 0;0 0 -1;0 1 0];
disp(s)
disp("The effect of multiplying S by a rigid body has the effect of swapping")
disp("the y and z indices and multiplying the z index by -1 as well as nullifying")
disp("the x index. There by aligning the rigid body to the Y-Z axes and flipping")
disp("the point pi radians (180-deg) about the X axis.")
disp("Example: [1 2 3] moves to [0 -3 2], moved inline with Y-Z axes,")
disp("flipped 180 about the X-axis.")
fprintf("\n")

%% 4
disp("----------------------------- Problems 4 --------------------------------")
syms phi_dot theta_dot psi_dot
phiDotxS = [0 0 0;0 0 -phi_dot;0 phi_dot 0];
disp("Multiplying S by the scalar 𝜙 dot yeilds")
disp(phiDotxS)
disp("This is the X component of the entire skew-symmetric angluar rate matrix ")
disp("W wrt the inertial frame.")
disp("For completeness Rdot_x is shown below:")
disp(phiDotxS*rotx)
fprintf("\n")

%% 5
disp("----------------------------- Problems 5 & 6 -----------------------------")
Rz = [cos(psi) -sin(psi) 0;sin(psi) cos(psi) 0;0 0 1];
Ry = [cos(theta) 0 sin(theta);0 1 0;-sin(theta) 0 cos(theta)];
Rx = [1 0 0;0 cos(phi) -sin(phi);0 sin(phi) cos(phi)];

ddtRz = [-sin(psi)*psi_dot -cos(psi)*psi_dot 0;cos(psi)*psi_dot -sin(psi)*psi_dot 0;0 0 0];
ddtRy = [-sin(theta)*theta_dot 0 cos(theta)*theta_dot;0 0 0;-cos(theta)*theta_dot 0 -sin(theta)*theta_dot];
ddtRx = [0 0 0; 0 -sin(phi)*phi_dot -cos(phi)*phi_dot; 0 cos(phi)*phi_dot -sin(phi)*phi_dot];

RzdotRyRx=ddtRz*Ry*Rx;
RydotRzRx=Rz*ddtRy*Rx;
RxdotRyRz=Rz*Ry*ddtRx;

disp("5. We are looking for some function that satisfies: Rdot = f(.)R so")
disp("that we may write w = T*edot.")
disp("first step: w = T*edot ---> Rdot = (w)^x * R, ")
disp("where the (.)^x is the skew systemric representation")
disp("Expanding Rdot: Rdot = Rzdot*Rz'+Rydot*Ry'+Rxdot*Rx'")
disp("From earlier we determined that Rxdot = phi_dot*S*Rx, I will refer to phi_dot*S as phi_dot^x")
disp("Now rearragning: (w)^x*R = ((psi_dot^x)Rz)RyRx + Rz((theta_dot^x)Ry)Rx + RzRy((phi_dot^x)Rx) ")
disp("Grouping and distributing: (w)^x*R = (psi_dot^x)R + Rz(theta_dot^x)R + RzRy(phi_dot^x)R ")
disp("Simplifying...")
disp("(w)^x*R = {(psi_dot^x)+Rz(theta_dot^x)+RzRy(phi_dot^x)}R")
disp("R*R(inv)=R*R'; simplifiying right-hand side...")

psi_dot_k = [0; 0 ;psi_dot];
theta_dot_j = [0; theta_dot; 0];
phi_dot_i = [phi_dot; 0; 0];

omega1=psi_dot_k;
omega2=Rz*theta_dot_j;
omega3=Rz*Ry*phi_dot_i;
omega=[omega1 omega2 omega3];
disp(omega)

disp("Extracting T")
T =[0 sin(psi) cos(psi)*cos(theta);0 cos(psi) cos(theta)*sin(psi);1 0 -sin(theta)];
disp(T)
disp("6. T is not a rotation matrix since no change in orientation to a frame occurs")
disp("after multiplying by T, simply a change of reference. T is idempotent on any frame A.")
fprintf("\n")

disp("------------------- Problems 7 -------------------")
fprintf("\n")
disp("The mechanical systems that I am considering for my project in this class are:")
disp("1. Virtual autonomous vehicle (2D representation)")
disp("2. Thrust Vector Controlled Rocket Motor/Engine")
disp("3. An Inverted pendulum (b/c I'm possibly swamped with other stuff)")
