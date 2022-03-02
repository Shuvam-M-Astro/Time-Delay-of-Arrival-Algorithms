import numpy as np
import matplotlib.pyplot as plt
import math
import sys
np.set_printoptions(threshold=sys.maxsize)
from numpy import loadtxt
#import copy
import warnings
from itertools import permutations
from itertools import combinations
from datetime import datetime
import time
from numpy import matlib as mb
from isqrt import isqrt
import random

from Fried import *
from Uncertainty_Fang import *
from Foy import *



warnings.filterwarnings("ignore", category=RuntimeWarning) 

############################################################################
############################################################################
############################################################################
############################################################################

R_set = 1e31


plt.close("all")

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


now = datetime.now()

current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

tic = time.perf_counter()  
#-----------------------------------------------------------------------------------

    
def Cart2Sph(x,y,z):
    if x == 0:
        x = 1e-100
    El =  np.arccos(z/(np.sqrt(x**2+y**2+z**2)))

    Az = np.arctan(y/x)
    R = np.sqrt(x**2+y**2+z**2)
    return R,El,Az

#Function to find x,y,z from spherical coordinates
def Sph2Cart(r,El,Az):
    x = r*np.cos(Az)*np.sin(El)
    y = r*np.sin(Az)*np.sin(El)
    z = r*np.cos(El)
    return x,y,z    

    
def Angsepcart(x1,y1,z1,x2,y2,z2):
    RA1 = Cart2Sph(x1,y1,z1)[2]
    Dec1 = np.pi/2- Cart2Sph(x1,y1,z1)[1]
    
    RA2 = Cart2Sph(x2,y2,z2)[2]
    Dec2 = Cart2Sph(x2,y2,z2)[1]
    
    Angsep1 = np.arccos(np.sin(Dec1)*np.sin(Dec2) + np.cos(Dec1)*np.cos(Dec2)*np.cos(RA1-RA2)           ) *rad
    return Angsep1#,Angsep2

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


# Checks if a matrix is a valid rotation matrix.
def isRotationMatrix(R) :
    Rt = np.transpose(R)
    shouldBeIdentity = np.dot(Rt, R)
    I = np.identity(3, dtype = R.dtype)
    n = np.linalg.norm(I - shouldBeIdentity)
    return n < 1e-6


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

def Az(Az):
        if Az > np.pi*2 :
            Az = Az - 2*np.pi
        elif   Az < 0  :
            Az = 2*np.pi + Az
            Az = Az
        return Az
        
def El(El):
    if El > np.pi/2:
        El = np.pi - El
    elif El < - np.pi/2:
        El = -np.pi - El
    return El

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
z1 = 1e-32

x2 = -0.866*r_E
y2 = -0.5*r_E
z2 = 1e-24

x3 = 1e-100
y3 = r_E
z3 = 1e-100

x4 = 1e-4
y4 = 1e-100
z4 = -r_E

x5 = 1e-100
y5 = -1*r_E
z5 = 1e-3

x6 = 0.866*r_E
y6 = 0.5*r_E
z6 = 1e-54

x7 = 0.866*r_E
y7 = -0.5*r_E
z7 = 1e-12

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

P_TPSP_LVLH_sel_0 = [ [0.127*r_E , 0.048*r_E , -0.041*r_E ,  0.030*r_E, -0.052*r_E, -0.136*r_E], 
 [0.272*r_E ,0.295*r_E, 0.299*r_E ,-0.296*r_E, -0.292*r_E, -0.267*r_E],[
 1e-32,  1e-3,  1e-21 , 1e-12 , 1e-34 , 1e-54] ]


P_TPSP_LVLH_opt_0 = [ [0.127*r_E , 0.048*r_E , -0.041*r_E  , 0.184*r_E, -0.132*r_E, -0.211*r_E],[ 
 0.272*r_E, 0.295*r_E, 0.299*r_E ,-0.237*r_E ,-0.269*r_E, -0.214*r_E],[
 1e-24 , 1e-6 , 1e-4 , 1e-23  ,1e-34 , 1e-54] ]

P_TPSP_LVLH_sel = []
P_TPSP_LVLH_opt = []

for i in range(0,len(P_TPSP_LVLH_sel_0[0])):
    P_TPSP_LVLH_sel.append([P_TPSP_LVLH_sel_0[0][i], P_TPSP_LVLH_sel_0[1][i],P_TPSP_LVLH_sel_0[2][i]    ]     )
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

############################################################################
############################################################################
############################################################################
############################################################################

def Err(val):

    RA = [val]# np.linspace(0,360,1000)
    
    val = 500
    DEC = np.hstack(  (     (np.linspace(91,179,val) ,   (np.linspace(1,89,val))              )         )                                )
    
    arr = DEC
    
    Err_Fried = []
    Err_Fried_2 = []
    Err_Fang = []
    Err_Foy_all = []
    Uncert_TPSP = [0]
    
    ########################### FRIEDLANDER ###########################################
    
    for i in range(0,len(RA)):
        for j in range(0,len(DEC)):
            P_TPSP = [D1,D2,D3,D5,D6,D7]
            Err_Fried.append(Friedlander_Err(P_TPSP,[0,0,0],RA[i],DEC[j]))
            
    ########################### FRIEDLANDER ###########################################
    
    ########################### FOY ###########################################
    
    P_TPSP = []
    P_TPSP.append(D1)
    P_TPSP.append(D2)
    P_TPSP.append(D3)
    P_TPSP.append(D5)
    P_TPSP.append(D6)
    P_TPSP.append(D7)
    
    for i in range(0,len(RA)):
        for j in range(0,len(DEC)):       
            Err_Foy_all.append(Foy(P_TPSP,Uncert_TPSP[0],RA[i],DEC[j],0,0)[0])
           
    ########################### FOY ###########################################
    
    
    ########################### FANG ###########################################
    
    P_TPSP = []
    P_TPSP.append(D1)
    P_TPSP.append(D2)
    P_TPSP.append(D3)
    P_TPSP.append(D5)
    P_TPSP.append(D6)
    P_TPSP.append(D7)

    Err_Fang = []
    
    for i in range(0,len(RA)):
        for j in range(0,len(DEC)):       
            if DEC[j] > 90:
                DEC[j] = DEC[j] - 180
            Fang_1 = Fang_Err(P_TPSP,[0,0,0],RA[i],DEC[j])
            Err_Fang.append(Fang_1)
           
    ########################### FANG ###########################################
    
    # Original DEC array
    DEC = np.hstack(  (     (np.linspace(91,179,val) ,   (np.linspace(1,89,val))              )         )                                )
    arr = DEC
    
    if arr[-1] == DEC[-1]:
        for i in range(0,len(arr)):
            if arr[i]>90:
                arr[i] = arr[i] - 180
                
                
                
    plt.figure()
    plt.plot(arr,Err_Foy_all,label='Foy_s Method with RA = ' +str(RA[0])+ ' degrees for TP/SP ')
    plt.plot(arr,Err_Fried,label='Friedlander_s Method with RA = ' +str(RA[0])+ ' degrees for TP/SP ')
    plt.plot(arr,Err_Fang,label='Fang_s Method with RA = ' +str(RA[0]) + ' degrees for TP/SP ')
    plt.legend(loc="upper left",fontsize='15')
    plt.yscale("log")
    plt.xlabel("DEC (Degrees)")
    plt.ylabel("abs(Err) (Degrees)")





Err_RA1 = Err(15)



















