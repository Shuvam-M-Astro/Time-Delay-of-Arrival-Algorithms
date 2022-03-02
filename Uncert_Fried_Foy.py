# Import Required Packages

import numpy as np
import matplotlib.pyplot as plt
import math
import sys
np.set_printoptions(threshold=sys.maxsize)
import warnings
from itertools import permutations
from itertools import combinations
from datetime import datetime
import time
from numpy import matlib as mb
from isqrt import isqrt
import random

from Foy import *
from Fried import *
# Ignore runtime warnings
warnings.filterwarnings("ignore", category=RuntimeWarning) 

############################################################################
############################################################################
############################################################################
############################################################################

# Close previous plots
plt.close("all")

# Set distance to the source
R_set = 1e31

#Define earth radius
r_E = 6.371e6 #m

#Convert radian to degrees
rad = 57.29577951308232

#Define the speed of light
c = 299792458;


############################################################################
############################################################################
############################################################################
############################################################################

#-----------------------------------------------------------------------------------
# Start timer
now = datetime.now()

current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

tic = time.perf_counter()  
#-----------------------------------------------------------------------------------

# Function to convert from Cartesian to Spherical coordinates    
def Cart2Sph(x,y,z):
    if x == 0:
        x = 1e-100
    DEC =  np.arccos(z/(np.sqrt(x**2+y**2+z**2)))

    RA = np.arctan(y/x)
    R = np.sqrt(x**2+y**2+z**2)
    return R,DEC,RA

# Function to convert from Spherical to Cartesian coordinates
def Sph2Cart(r,DEC,RA):
    x = r*np.cos(RA)*np.sin(DEC)
    y = r*np.sin(RA)*np.sin(DEC)
    z = r*np.cos(DEC)
    return x,y,z    

# Function to find the angular separation of two points using their cartesian coordinates
def Angsepcart(x1,y1,z1,x2,y2,z2):
    RA1 = Cart2Sph(x1,y1,z1)[2]
    Dec1 = np.pi/2- Cart2Sph(x1,y1,z1)[1]
    
    RA2 = Cart2Sph(x2,y2,z2)[2]
    Dec2 = Cart2Sph(x2,y2,z2)[1]
    
    Angsep1 = np.arccos(np.sin(Dec1)*np.sin(Dec2) + np.cos(Dec1)*np.cos(Dec2)*np.cos(RA1-RA2)           ) *rad
    return Angsep1#,Angsep2

# Function to find the angular separation of two points using their spherical coordinates
def Angsepsph(RA1,Dec1,RA2,Dec2):
    Angsep1 = np.arccos(np.sin(Dec1)*np.sin(Dec2) + np.cos(Dec1)*np.cos(Dec2)*np.cos(RA1-RA2)           ) *rad
    return Angsep1

#Function to find rotation matrix for a given set of equations
def rot(vec1, vec2):
   
    a, b = (vec1 / np.linalg.norm(vec1)).reshape(3), (vec2 / np.linalg.norm(vec2)).reshape(3)
    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)
    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    rotation_matrix = np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))

    return rotation_matrix


# Function to checks if a matrix is a valid rotation matrix.
def isRotationMatrix(R) :
    Rt = np.transpose(R)
    shouldBeIdentity = np.dot(Rt, R)
    I = np.identity(3, dtype = R.dtype)
    n = np.linalg.norm(I - shouldBeIdentity)
    return n < 1e-6

# Function to find Euler Angles from rotation matrix
def rotationMatrixToEulerAngles(R) :
    assert(isRotationMatrix(R))
    sy = math.sqrt(R[0,0] * R[0,0] +  R[1,0] * R[1,0])
    singular = sy < 1e-6

    if  not singular :
        x = math.atan2(R[2,1] , R[2,2])
        y = math.atan2(-R[2,0], sy)
        z = math.atan2(R[1,0], R[0,0])

    else :
        x = math.atan2(-R[1,2], R[1,1])
        y = math.atan2(-R[2,0], sy)
        z = 0

    return np.array([x, y, z])

# Constrain RA between 0 - 2*pi
def RA(RA):
        if RA > np.pi*2 :
            RA = RA - 2*np.pi
        elif   RA < 0  :
            RA = 2*np.pi + RA
            RA = RA
        return RA
 
# Constrain DEC between 0 - np.pi/2       
def DEC(DEC):
    if DEC > np.pi/2:
        DEC = np.pi - DEC
    elif DEC < - np.pi/2:
        DEC = -np.pi - DEC
    return DEC

def length(x1,y1,z1,x2,y2,z2):
   A = np.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
   return A

def Angsep2(RA1,RA2,Dec1,Dec2):
    Angsep1 = np.arccos(np.sin(Dec1)*np.sin(Dec2) + np.cos(Dec1)*np.cos(Dec2)*np.cos(RA1-RA2)           ) *rad
    #Angsep2 = np.sqrt( ((RA1-RA2)*np.cos(Dec1))**2 + (Dec1-Dec2)**2            )*rad
    return Angsep1#,Angsep2


def more(N):
    Det = []
    
    A_all_1000 = []
    B_all_1000 = []
    C_all_1000 = []
    
    
    for i in range(0,N):
        A = random.randint(0,6.371e6)
        A = random.choice([-A, A])
        A_left = np.sqrt( r_E**2 - A**2           )
        B = random.randint(0,int(A_left))
        B = random.choice([-B, B])
        B_left = np.sqrt( r_E**2 - A**2   - B**2         )
        C = -B_left
        C = random.choice([-B_left, B_left])
        #C = random.randint(-int(B_left),int(B_left))
        
        A_all_1000.append(A)
        B_all_1000.append(B)
        C_all_1000.append(int(C))
    
        
        Det.append([A,B,C])
    
    
    return Det


def moreEq(N):
    Det = []
    
    A_all_1000 = []
    B_all_1000 = []
    C_all_1000 = []
    
    
    for i in range(0,N):
        A = random.randint(0,6.371e6)
        A = random.choice([-A, A])
        A_left = np.sqrt( r_E**2 - A**2           )
        
        B = random.choice([-A_left, A_left])
        #B = random.choice([-B, B])
        B_left = np.sqrt( r_E**2 - A**2   - B**2         )
        
        
        A_all_1000.append(A)
        B_all_1000.append(B)
        C = 1e-100
        Det.append([A,B,C])
    
    
    return Det


####################################################################
####################################################################
######################     Loop         ############################
####################################################################
####################################################################

#Detector configuration

x1 = -0.866*r_E
y1 = 0.5*r_E # + 1e3
z1 = 1e-6

x2 = -0.866*r_E
y2 = -0.5*r_E
z2 = 1e-80

x3 = 1e-100
y3 = r_E
z3 = 1e-95

x4 = 1e-4
y4 = 1e-100
z4 = -r_E

x5 = 1e-100
y5 = -1*r_E
z5 = 1e-90

x6 = 0.866*r_E
y6 = 0.5*r_E
z6 = 1e-85

x7 = 0.866*r_E
y7 = -0.5*r_E
z7 = 1e-100

x8 = 1e-2
y8 = 1e-100
z8 = r_E

x9 = r_E*np.cos(50/rad) 
z9 = np.sqrt(r_E**2-x9**2)
y9 = 1e-2

x10 = -r_E*np.sin(40/rad)
z10 = np.sqrt(r_E**2-x10**2)
y10 = 1e-3

z11 = -r_E*np.cos(40/rad)
x11 = -1*np.sqrt(r_E**2-z11**2)
y11= 2e-3

z12 = -r_E*np.sin(40/rad)
x12 = np.sqrt(r_E**2-z12**2)
y12 = 3e2

# 156 for eq triangle
# 123 for one sided detectors

P_TPSP_LVLH_sDEC_0 = [ [0.127*r_E , 0.048*r_E , -0.041*r_E ,  0.030*r_E, -0.052*r_E, -0.136*r_E], 
 [0.272*r_E ,0.295*r_E, 0.299*r_E ,-0.296*r_E, -0.292*r_E, -0.267*r_E],[
 1e-32,  1e-3,  1e-21 , 1e-12 , 1e-34 , 1e-54] ]


P_TPSP_LVLH_opt_0 = [ [0.127*r_E , 0.048*r_E , -0.041*r_E  , 0.184*r_E, -0.132*r_E, -0.211*r_E],[ 
 0.272*r_E, 0.295*r_E, 0.299*r_E ,-0.237*r_E ,-0.269*r_E, -0.214*r_E],[
 1e-24 , 1e-6 , 1e-4 , 1e-23  ,1e-34 , 1e-54] ]


P_TPSP_LVLH_sDEC = []
P_TPSP_LVLH_opt = []



for i in range(0,len(P_TPSP_LVLH_sDEC_0[0])):
    P_TPSP_LVLH_sDEC.append([P_TPSP_LVLH_sDEC_0[0][i], P_TPSP_LVLH_sDEC_0[1][i],P_TPSP_LVLH_sDEC_0[2][i]    ]     )
    P_TPSP_LVLH_opt.append([P_TPSP_LVLH_opt_0[0][i], P_TPSP_LVLH_opt_0[1][i],P_TPSP_LVLH_opt_0[2][i]    ]     )


X1 = [x1,y1,z1]
X2 = [x2,y2,z2]
X3 = [x3,y3,z3]
X4 = [x4,y4,z4]
X5 = [x5,y5,z5]
X6 = [x6,y6,z6]
X7 = [x7,y7,z7]
X8 = [x8,y8,z8]
X9 = [x9,y9,z9]
X10 = [x10,y10,z10]
X11 = [x11,y11,z11]
X12 = [x12,y12,z12]

D1 = X1
D2 = X2
D3 = X3
D4 = X4
D5 = X5
D6 = X6
D7 = X7
D8 = X8
D9 = X9
D10 = X10
D11 = X11
D12 = X12


def Uncertainty_all(Source,P,Width):
    
    P1 = np.copy(P)
    Fin_RA = []
    Fin_DEC = []
    
    RA = Source[0][0]
    DEC = Source[0][1]
    
    Fried1 = Uncertainty_Friedlander(RA,DEC,P1,Width)   
    Foy1 = Uncertainty_Foy(RA,DEC,P1,Width)   
    
    max_Err = min(abs(Fried1),abs(Foy1))
    
    Err = np.arange(0,max_Err,max_Err/50)
    
    RA_min,RA_max = RA-max_Err,RA+max_Err
    DEC_min,DEC_max = DEC-max_Err,DEC+max_Err, 
   
    for r in range(0,len(Err)):
        for a in range(0,360):
            Fin_RA.append(RA + Err[r] * np.cos(a))
            Fin_DEC.append(DEC + Err[r] * np.sin(a))
              
    return Fin_RA,Fin_DEC

##############################################################################
##############################################################################
######################         Uncertainty PLot          #####################
##############################################################################
##############################################################################

P = []

P.append(D1)
P.append(D2)
P.append(D3)

P.append(D5)
P.append(D6)
P.append(D7)


Fin_RA = []
Fin_DEC = []

RA = [165,15,115,65]
DEC = [15,45,75,-15,-45,-75]

Source = [[[15,15]],[[15,45]],[[15,75]],[[15,-15]],[[15,-45]],[[15,-75]],[[65,15]],[[65,45]],[[65,75]],[[65,-15]],[[65,-45]],[[65,-75]],[[115,15]],[[115,45]],[[115,75]],[[115,-15]],[[115,-45]],[[115,-75]],[[165,15]],[[165,45]],[[165,75]],[[165,-15]],[[165,-45]],[[165,-75]],              ]


for i in range(0,len(Source)):
    Uncert_All = Uncertainty_all(Source[i],P,1e-3)   
    Fin_RA.append(np.divide(Uncert_All[0],rad))
    Fin_DEC.append(np.divide(Uncert_All[1],rad))

Fin_RA_F = []
Fin_DEC_F = [] 

for i in range(0,len(Fin_RA)):
    for j in range(0,len(Fin_RA[i])):
        Fin_RA_F.append(Fin_RA[i][j])

for i in range(0,len(Fin_DEC)):
    for j in range(0,len(Fin_DEC[i])):
        Fin_DEC_F.append(Fin_DEC[i][j])

######################################################


Fin_RA2 = []
Fin_DEC2 = []

for i in range(0,len(Source)):
    Uncert_All = Uncertainty_all(Source[i],P,1e-4)   
    Fin_RA2.append(np.divide(Uncert_All[0],rad))
    Fin_DEC2.append(np.divide(Uncert_All[1],rad))

Fin_RA_F2 = []
Fin_DEC_F2 = [] 

for i in range(0,len(Fin_RA2)):
    for j in range(0,len(Fin_RA2[i])):
        Fin_RA_F2.append(Fin_RA2[i][j])

for i in range(0,len(Fin_DEC2)):
    for j in range(0,len(Fin_DEC2[i])):
        Fin_DEC_F2.append(Fin_DEC2[i][j])
     
fig = plt.figure(1,dpi=300)
ax = fig.add_subplot(111, projection="hammer")
#• here 111 means “subplot 1 of a 1x1 grid of plots”

#Set marker size
size = 1

R  = 1
x = []
y = []
z = []

ax.scatter(Fin_RA_F, Fin_DEC_F,marker='o',s=1,alpha=0.1,label='Δδt_ij = 1e-3s')
ax.scatter(Fin_RA_F2, Fin_DEC_F2,marker='o',s=1,alpha=0.1,label='Δδt_ij = 1e-4s')

RA_1 = []
DEC_1 = []
True_S = []

for i in range(0,len(RA)):
    for j in range(0,len(DEC)):
        True_S.append([RA[i],DEC[j]])
for i in range(0,len(True_S)):
    RA_1.append(True_S[i][0])
    DEC_1.append(True_S[i][1])   
 
ax.scatter(np.divide(RA_1,rad),np.divide(DEC_1,rad),marker='x',color='r',s=0.5,alpha=1,label='True position');  

ax.legend(loc='upper left')
xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
ax.set_xticklabels(xlab, weight=8,fontsize=4)
plt.legend(loc='upper left',fontsize=4,facecolor='white', framealpha=1)
ax.grid(color='k', linestyle='dotted', linewidth=0.5)

#################### TIME THE SCRIPT #############################
toc = time.perf_counter()  
print("Total time = " +str(toc - tic) +"seconds")
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)



























