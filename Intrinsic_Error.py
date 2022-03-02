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
    DEC =  np.arccos(z/(np.sqrt(x**2+y**2+z**2)))

    RA = np.arctan(y/x)
    R = np.sqrt(x**2+y**2+z**2)
    return R,DEC,RA

#Function to find x,y,z from spherical coordinates
def Sph2Cart(r,DEC,RA):
    x = r*np.cos(RA)*np.sin(DEC)
    y = r*np.sin(RA)*np.sin(DEC)
    z = r*np.cos(DEC)
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

def RA(RA):
        if RA > np.pi*2 :
            RA = RA - 2*np.pi
        elif   RA < 0  :
            RA = 2*np.pi + RA
            RA = RA
        return RA
        
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
z8 = -r_E

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


P_TPSP_LVLH_sDEC = [ [0.127*r_E , 0.048*r_E , -0.041*r_E ,  0.030*r_E, -0.052*r_E, -0.136*r_E], 
 [0.272*r_E ,0.295*r_E, 0.299*r_E ,-0.296*r_E, -0.292*r_E, -0.267*r_E],[
 1e-100,  1e-100,  1e-100 , 1e-100 , 1e-100 , 1e-100] ]


P_TPSP_LVLH_opt = [ [0.127*r_E , 0.048*r_E , -0.041*r_E  , 0.184*r_E, -0.132*r_E, -0.211*r_E],[ 
 0.272*r_E, 0.295*r_E, 0.299*r_E ,-0.237*r_E ,-0.269*r_E, -0.214*r_E],[
 1e-100 , 1e-100 , 1e-100 , 1e-100  ,1e-100 , 1e-100] ]




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

def Friedlander(P,Err,RA_set,DEC_set):   

    #Function for constraints on RAimuth and altitude

        def RA(RA):
            if RA > np.pi*2 :
                RA = RA - 2*np.pi
            elif   RA < 0  :
                RA = 2*np.pi + RA
                RA = RA
            return RA
        
        def DEC(DEC):
            if DEC > np.pi/2:
                DEC = np.pi - DEC
            elif DEC < - np.pi/2:
                DEC = -np.pi - DEC
            return DEC
    
    
        RA_set = RA_set/rad
        DEC_set = DEC_set/rad
          
        D1 = P[0]
        D2 = P[1]
        D3 = P[2]
        D4 = P[3]
        
        x1 = D1[0]
        y1 = D1[1]
        z1 = D1[2]
        
        x2 = D2[0]
        y2 = D2[1]
        z2 = D2[2]
        
        x3 = D3[0]
        y3 = D3[1]
        z3 = D3[2]
        
        x4 = D4[0]
        y4 = D4[1]
        z4 = D4[2]
 
        x1_1 = x1
        y1_1 = y1
        z1_1 = z1
        
        # Translate detector configuration
        
        x1 = x1 - x1_1 + 1e-100
        y1 = y1 - y1_1  
        z1 = z1 - z1_1
   
        x2 = x2 - x1_1 
        y2 = y2 - y1_1
        z2 = z2 - z1_1

        x3 = x3 - x1_1
        y3 = y3 - y1_1
        z3 = z3 - z1_1

        x4 = x4 - x1_1
        y4 = y4 - y1_1
        z4 = z4 - z1_1

        R_set = 1e21
        Source = Sph2Cart(R_set, DEC_set, RA_set)
        xs = Source[0]
        ys = Source[1]
        zs = Source[2]
 
        xs = xs - x1_1
        ys = ys - y1_1
        zs = zs - z1_1
        
        #Define the distances between the detectors
        d21 = np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
        d31 = np.sqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)
        d41 = np.sqrt((x4-x1)**2+(y4-y1)**2+(z4-z1)**2)
        
        d12 = d21
        d13 = d31
        d14 = d41
        
        #Calculate Spherical coordinates from cartesian coordinates ... Only for angular distance
        
        RAs = RA_set
        DECs = DEC_set
        
        RAar = RAs
        DECar = np.pi/2 - DECs
        DECar = DEC(DECar)
        RAar = RA(RAar)
        
        R1  = Cart2Sph(x1,y1,z1)[0]
        DEC1 = np.pi/2 - Cart2Sph(x1,y1,z1)[1]
        RA1 = Cart2Sph(x1,y1,z1)[2]
        DEC1 = DEC(DEC1)
        RA1 = RA(RA1)
        
        R2  = Cart2Sph(x2,y2,z2)[0]
        DEC2 = np.pi/2 - Cart2Sph(x2,y2,z2)[1]
        RA2 = Cart2Sph(x2,y2,z2)[2]
        DEC2 = DEC(DEC2) 
        RA2 = RA(RA2)
        
        R3  = Cart2Sph(x3,y3,z3)[0]
        DEC3 = np.pi/2 - Cart2Sph(x3,y3,z3)[1]
        RA3 = Cart2Sph(x3,y3,z3)[2]
        DEC3 = DEC(DEC3)
        RA3 = RA(RA3)
        
        R4  =           Cart2Sph(x4,y4,z4)[0]
        DEC4 = np.pi/2 - Cart2Sph(x4,y4,z4)[1]
        RA4 =           Cart2Sph(x4,y4,z4)[2]
        DEC4 = DEC(DEC4)
        RA4 = RA(RA4)
        
        # Calculate time-dDECays from RA and DEC
        
        [xar,yar,zar] = [xs,ys,zs]
        
        d1ar = isqrt( (int(float(format(xar, '.100g'))) - int(float(format(x1, '.100g'))))**2                   + (int(float(format(yar, '.100g')) - int(float(format(y1, '.100g'))))**2       + (int(float(format(zar, '.100g')))          - int(float(format(z1, '.100g')))))**2  )  
               
        d2ar = isqrt( (int(float(format(xar, '.100g'))) - int(float(format(x2, '.100g'))))**2                   + (int(float(format(yar, '.100g')) - int(float(format(y2, '.100g'))))**2       + (int(float(format(zar, '.100g')))          - int(float(format(z2, '.100g')))))**2  ) 
     
        d3ar = isqrt( (int(float(format(xar, '.100g'))) - int(float(format(x3, '.100g'))))**2                   + (int(float(format(yar, '.100g')) - int(float(format(y3, '.100g'))))**2       + (int(float(format(zar, '.100g')))          - int(float(format(z3, '.100g')))))**2  ) 

        d4ar = isqrt( (int(float(format(xar, '.100g'))) - int(float(format(x4, '.100g'))))**2                   + (int(float(format(yar, '.100g')) - int(float(format(y4, '.100g'))))**2       + (int(float(format(zar, '.100g')))          - int(float(format(z4, '.100g')))))**2  ) 
        
        Z_all = []
        for i in range(0,len(P)):
            Z_all.append(P[i][2])
        
        t21ar = (d2ar-d1ar)/c 
        t21 = t21ar; 
        t12 = -t21 + Err[0]
        
        t31ar = (d3ar-d1ar)/c 
        t31 = t31ar;   
        t13 = -t31 + Err[1]
    
        t41ar = (d4ar-d1ar)/c 
        t41 = t41ar;   
        t14 = -t41 + Err[2]
            
        
        Ang2s = Angsepsph(RA2,DEC2,RAar,DECar)/rad
        Ang3s = Angsepsph(RA3,DEC3,RAar,DECar)/rad
        Ang4s = Angsepsph(RA4,DEC4,RAar,DECar)/rad
         
        
        # Time-delay and add uncertainty
        t12 = np.cos(Ang2s) * d12 / c 
   
        t12 = t12 + Err[0]
         
        t13 = np.cos(Ang3s) * d13 / c 
        t13 = t13 + Err[1]
         
        t14 = np.cos(Ang4s) * d14 / c 
        t14 = t14 + Err[2]
    
        # Left side and right side
    
        S_12_Z = np.array([x1-x2,y1-y2,z1-z2]).transpose()
        S_13_Z = np.array([x1-x3,y1-y3,z1-z3]).transpose()
        S_14_Z = np.array([x1-x4,y1-y4,z1-z4]).transpose()
        
        L = np.array([S_12_Z  , S_13_Z,  S_14_Z])
       
        R =  np.array([[d12*t12/(d12/c)],[d13*t13/(d13/c)],[d14*t14/(d14/c)]]  )
        
        U = np.linalg.solve(L,R)
     
        U[0] = U[0] 
        U[1] = U[1] 
        U[2] = U[2] 
        
        RA = Cart2Sph(U[0],U[1],U[2])[2]*rad
        for item in RA:
            RA = float(RA)
            
        if RA<0:
            RA = RA + np.pi*rad
        
        DEC_calc = Cart2Sph(U[0],U[1],U[2])[1]*rad
        
        if RAar > 180/rad:
            RA = 180+RA
        
        if DEC_set>0:
            DEC_calc = np.pi*rad - Cart2Sph(U[0],U[1],U[2])[1]*rad
            for item in DEC_calc:
                DEC_calc = float(DEC_calc)
        
        elif DEC_set<0:       
            DEC_calc =  Cart2Sph(U[0],U[1],U[2])[1]*rad - np.pi*rad
            
            for item in DEC_calc:
                DEC_calc = float(DEC_calc)
        
        Tot_er = abs(Angsepsph(RA_set,DEC_set,RA/rad,DEC_calc/rad))
        
        return Tot_er,RA_set*rad,DEC_set*rad,RA,DEC_calc

    
def FangD3(D1,D2,D3,Err,RA_set,DEC_set):
    if DEC_set<90:
        DEC_set = 90 - DEC_set
    
    #################################################
    x1 = D1[0]
    y1 = D1[1]
    z1 = D1[2]
    
    x2 = D2[0]
    y2 = D2[1]
    z2 = D2[2]

    x3 = D3[0]
    y3 = D3[1]
    z3 = D3[2]

    x1_1 = x1
    y1_1 = y1
    z1_1 = z1
    
    # Translate detector configuration
    
    x1 = x1 - x1_1 + 1e-6
    y1 = y1 - y1_1
    z1 = z1 - z1_1

    x2 = x2 - x1_1 
    y2 = y2 - y1_1
    z2 = z2 - z1_1
    
    x3 = x3 - x1_1
    y3 = y3 - y1_1
    z3 = z3 - z1_1

    #Assign coordinate matrices
    X1 = [x1,y1,z1]
    X2 = [x2,y2,z2]
    X3 = [x3,y3,z3]
    #Define a source position in spherical coordinates and find the cartesian coordinates
    
    #############################################

    P = [D1,D2,D3]
    
    P_first = [P[0][0],P[0][1],P[0][2]]
    for i in range(0,len(P)):
            P[i][0],P[i][1],P[i][2] = P[i][0] - P_first[0]  , P[i][1] - P_first[1], P[i][2]  - P_first[2]

    DECs =  DEC_set  #Top Down (0 - 180)
    RAs =  RA_set/rad # 0 - 360

    #Calculate Spherical coordinates from cartesian coordinates
    
    DEC_all = []
    RA_all = []
    
    for i in range(0,len(P)):
        DEC_all.append(DEC( np.pi/2 -Cart2Sph(P[i][0],P[i][1],P[i][2])[1]   ))
        RA_all.append(RA(Cart2Sph(P[i][0],P[i][1],P[i][2])[2]   ))
   
    DECar = np.pi/2 - DECs/rad
    RAar = RAs/rad
    DECar = DEC(DECar)
    RAar = RA(RAs)
    
    
    ######################################################
    #R1 = Cart2Sph(x1,y1,z1)[0]
    DEC1 = np.pi/2 - Cart2Sph(x1,y1,z1)[1]
    RA1 = Cart2Sph(x1,y1,z1)[2]
    DEC1 = DEC(DEC1)
    RA1 = RA(RA1)

    #R2 = Cart2Sph(x2,y2,z2)[0]
    DEC2 = np.pi/2 - Cart2Sph(x2,y2,z2)[1]
    RA2 = Cart2Sph(x2,y2,z2)[2]
    DEC2 = DEC(DEC2) 
    RA2 = RA(RA2)

    #R3 = Cart2Sph(x3,y3,z3)[0]
    DEC3 = np.pi/2 - Cart2Sph(x3,y3,z3)[1]
    RA3 = Cart2Sph(x3,y3,z3)[2]
    DEC3 = DEC(DEC3)
    RA3 = RA(RA3)
    ###############################################################

    ########## Spherical coordinates of origin
    DEC0 = 0
    RA0 = 0
      
    #-------------------------------------------------------------------------------------
    
    #Define the distances between the detectors
    d_all = [] 
    for i in range(0,len(P)-1):
            d_all.append(length(P[0][0],P[0][1],P[0][2],P[i+1][0],P[i+1][1],P[i+1][2]))
    
    ###########################################
    
    d21 = np.sqrt((x2-x1)**2+(y2-y1)**2+(z2-z1)**2)
    d31 = np.sqrt((x3-x1)**2+(y3-y1)**2+(z3-z1)**2)
    d12=d21
    d13=d31
   
###########################################
    # Distances between the detectors
    r_all = []
    for i in range(0,len(P)-1):
        r_all.append( [P[i+1][0] - P[0][0],P[i+1][1] - P[0][1],P[i+1][2] - P[0][2]     ]                         )

    ##############################################################################################
    ###############################      TRANSFORMATION              #############################
    ##############################################################################################

    #######################################-Rotation #################################
    
    
    ######## Find the coordinates of the rotation matrix from the positions of the detectors

    r11 = r_all[0][0]
    r12 = r_all[0][1]
    r13 = r_all[0][2]

    r21 = r_all[1][0]
    r22 = r_all[1][1]
    r23 = r_all[1][2]
    
    ri = r12*r23-r22*r13
    rj = r11*r23-r21*r13 
    rk = r11*r22-r21*r12
    
    a = ri
    b = rj
    c = rk
    d = ri*x1+rj*y1+rk*z1
    
    #Find the rotation matrix using the normal vectors of the planes
    
    vector_1 = [a,b,c]
    vector_2 = [0,1,0]
    unit_vector_1 = vector_1 / np. linalg. norm(vector_1)
    unit_vector_2 = vector_2 / np. linalg. norm(vector_2)
    dot_product = np. dot(unit_vector_1, unit_vector_2)
    angle = np. arccos(dot_product)

    axis = np.cross(vector_1,vector_2)/np.linalg.norm(np.cross(vector_1,vector_2))

    cos_angle = np.cos(angle)
    s = np.sqrt(1-cos_angle*cos_angle)
    C = 1-cos_angle
    
    x = axis[0]
    y = axis[1]
    z = axis[2]
    
    #Matrix DECements for rotation matrix
    r_mat_11 = x*x*C+cos_angle
    r_mat_12 = x*y*C-z*s
    r_mat_13 = x*z*C+y*s
    
    r_mat_21 = y*x*C+z*s
    r_mat_22 = y*y*C+cos_angle
    r_mat_23 = y*z*C-x*s
    
    r_mat_31 = z*x*C-y*s
    r_mat_32 = z*y*C+x*s
    r_mat_33 = z*z*C+cos_angle
    
    #Final rotation matrix
    rmat = [[r_mat_11,r_mat_12,r_mat_13],[r_mat_21,r_mat_22,r_mat_23],[r_mat_31,r_mat_32,r_mat_33]]
    
    #Find detector coordinates in new frame
    
    #Rotation
    X1_2  = np.matmul(np.linalg.inv(rmat),X1) 
    X2_2  = np.matmul(np.linalg.inv(rmat),X2) 
    X3_2  = np.matmul(np.linalg.inv(rmat),X3)
    
    # Convert to spherical coordinates
    
    DEC1_2 = np.pi/2 - Cart2Sph(X1_2[0],X1_2[1],X1_2[2])[1]
    DEC1_2 = DEC(DEC1_2)
    RA1_2 = Cart2Sph(X1_2[0],X1_2[1],X1_2[2])[2]
    RA1_2 = RA(RA1_2)
    
    DEC2_2 = np.pi/2 - Cart2Sph(X2_2[0],X2_2[1],X2_2[2])[1]
    DEC2_2 = DEC(DEC2_2)
    RA2_2 = Cart2Sph(X2_2[0],X2_2[1],X2_2[2])[2]
    RA2_2 = RA(RA2_2)
    
    DEC3_2 = np.pi/2 - Cart2Sph(X3_2[0],X3_2[1],X3_2[2])[1]
    DEC3_2 = DEC(DEC3_2)
    RA3_2 = Cart2Sph(X3_2[0],X3_2[1],X3_2[2])[2]
    RA3_2 = RA(RA3_2)
    
    # Compare angles in the new frame to those in the original frame
    
    d12_2 = length(X1_2[0],X1_2[1],X1_2[2],X2_2[0],X2_2[1],X2_2[2])
    d13_2 = length(X1_2[0],X1_2[1],X1_2[2],X3_2[0],X3_2[1],X3_2[2])
    
    #Rotation matrices to direct euler angles
    
    rmat = np.array(np.linalg.inv(rmat))
    isRotationMatrix(rmat)
  
    if isRotationMatrix(rmat) == True:
    
        E_1 = rotationMatrixToEulerAngles(rmat)
        

        RAar_2 = RAar + E_1[1]
        DECar_2 = DECar + E_1[2]
        

        ##########################################################################
        ####################  Time-dDECays from RA and DEC ########################
        ##########################################################################
        R_set = 1e31
        [xs,ys,zs] = [Sph2Cart(R_set,DECar_2,RAar_2)[0],Sph2Cart(R_set,DECar_2,RAar_2)[1],Sph2Cart(R_set,DECar_2,RAar_2)[2]]
    
        [xar,yar,zar] = [xs,ys,zs]
    
        R_set = 1e6
    
        [x1,y1,z1] = [X1_2[0],X1_2[1],X1_2[2]]
        [x2,y2,z2] = [X2_2[0],X2_2[1],X2_2[2]]
        [x3,y3,z3] = [X3_2[0],X3_2[1],X3_2[2]]
        
        
        d1ar = (isqrt( (int(float(format(xar, '.100g'))) - int(float(format(x1, '.100g'))))**2                   + (int(float(format(yar, '.100g'))) - int(float(format(y1, '.100g'))))**2       + (int(float(format(zar, '.100g')))          - int(float(format(z1, '.100g'))))**2  ) )

        d2ar = (isqrt( (int(float(format(xar, '.100g'))) - int(float(format(x2, '.100g'))))**2                   + (int(float(format(yar, '.100g'))) - int(float(format(y2, '.100g'))))**2       + (int(float(format(zar, '.100g')))          - int(float(format(z2, '.100g'))))**2  ) )

        d3ar = (isqrt( (int(float(format(xar, '.100g'))) - int(float(format(x3, '.100g'))))**2                   + (int(float(format(yar, '.100g'))) - int(float(format(y3, '.100g'))))**2       + (int(float(format(zar, '.100g')))          - int(float(format(z3, '.100g'))))**2  ) )
        
      
        c = 299792458

        # Time-dDECay calculation and add error
        t21ar = (d2ar-d1ar)/c 
        t21 = t21ar; 
        t12 = -t21 + Err[0]
        
        t31ar = (d3ar-d1ar)/c 
        t31 = t31ar;   
        t13 = -t31

        # Angles between pairs of detectors w.r.t D1
        Br12 = np.arccos(299792458*t12/d12)*rad 
        Br13 = np.arccos(299792458*t13/d13)*rad 
        
        ########### Find the source coordinates in new frame
        DDECs_2 =  np.arcsin(   ((np.cos(Br12/rad) - (np.cos(DEC2_2) / np.cos(DEC3_2)) * np.cos(Br13/rad)) / (np.sin(DEC2_2) - (np.cos(DEC2_2) / np.cos(DEC3_2)) * np.sin(DEC3_2)))) 
        RAs_2  =  RA2_2 + np.arccos(((np.cos(Br12/rad) - np.sin(DEC2_2) * np.sin(DDECs_2)) / (np.cos(DEC2_2) * np.cos(DDECs_2))))
        RAs_2_2 = np.pi*2 - RAs_2
        
        #Errors in 2nd frame
        
        RA2_1_er = (abs((RAs_2*rad  - RAar_2*rad)/(RAar_2*rad))*100)
        RA2_2_er = (abs((np.pi*2*rad - RAs_2*rad  - RAar_2*rad)/(RAar_2*rad))*100)
        DEC2_er = (abs((DDECs_2*rad - DECar_2*rad)/(DECar_2*rad))*100)
        
        # Rotation to original Frame of reference
        FinRA_est =  180 -  RAs_2*rad - E_1[1]*rad  
        FinDEC_est = np.pi/2*rad -( 1 * (DDECs_2*rad - E_1[2]*rad)) 
            
        # Constrain calculated estimate
        if FinDEC_est < np.pi/2*rad:  
            FinDEC_est = FinDEC_est
        if FinDEC_est > np.pi/2*rad:
            FinDEC_est = np.pi/2*rad -( 1 * (DDECs_2*rad - E_1[2]*rad)) - np.pi*rad #- E_2[2]*rad)
        
        if DECar_2<0:
            DECar_2 = -1*DECar_2
        
        if RAar_2>180/rad:
            FinRA_est = 360 - FinRA_est
            
        if DEC_set<90:
            FinRA_est = 180 - FinRA_est

        Fin_Err = abs(Angsep2(FinRA_est/rad,RAar_2,FinDEC_est/rad,DECar_2))
        
        if DECar_2*rad<90:
            FinDEC_est = 90 - FinDEC_est
            
        if DECar_2>90/rad:
            FinDEC_est = 0

        if FinRA_est > 360:
            FinRA_est  = FinRA_est - 360
        
        if FinRA_est < 0:
            FinRA_est  = FinRA_est + 360*2
        
        
        return Fin_Err,FinRA_est,FinDEC_est,RAar_2*rad,DECar_2*rad
    else:
        # In case the rotation matrix was uninvertible
        return [1e6,1e6,1e6,1e6,1e6]

    








def Foy(Pini,Err,RAs,DECs):
        
        Iter = 35

        r = R_set
        
        P = np.copy(Pini)
        Err_all = [1000,999]
        RA_all = [1000,999]
        DEC_all = [1000,999]
        
        P_first = [P[0][0],P[0][1],P[0][2]]
        for i in range(0,len(P)):
            P[i][0],P[i][1],P[i][2] = P[i][0] - P_first[0]  , P[i][1] - P_first[1], P[i][2]  - P_first[2]

        M = len(P)
        
        RAs = RAs/rad 
        DECs = DECs /rad  
        
        xs = r*np.cos(RAs)*np.sin(DECs)
        ys = r*np.sin(RAs)*np.sin(DECs)
        zs = r*np.cos(DECs)
         
        xs = xs - P_first[0]
        ys = ys - P_first[1]
        zs = zs - P_first[2]
        
        p_T = [[xs], [ys ], [zs]]
        
        Ang = []
        for i in range(0,len(P)-1):
            Ang.append(Angsepcart(P[i+1][0],P[i+1][1],P[i+1][2],xs,ys,zs)      )

        d_all = [] 
        for i in range(0,len(P)-1):
            d_all.append(length(P[0][0],P[0][1],P[0][2],P[i+1][0],P[i+1][1],P[i+1][2]))
        
        toa_1 = []
        dummy = mb.repmat(p_T,1,M)    -    np.transpose(P)
        
        for ii in range(0,M):
           toa_1.append(np.linalg.norm(np.array(dummy[:,ii]))/c    )
           
        tdoa = []
        for i in range(1,len(toa_1)):
          tdoa.append(toa_1[i] - toa_1[0])
          
        d_all = []
        [xar,yar,zar] = [xs,ys,zs]
        
        for i in range(0,len(Pini)):
            d_all.append(isqrt( (int(float(format(xar, '.100g'))) - int(float(format(Pini[i][0], '.100g'))))**2                   + (int(float(format(yar, '.100g'))) - int(float(format(Pini[i][1], '.100g'))))**2       + (int(float(format(zar, '.100g')))          - int(float(format(Pini[i][2], '.100g'))))**2  ) )
         
        t_all = []
        
        for i in range(1,len(d_all)):
            t_all.append( (d_all[i] - d_all[0]      ) /c            )
        
        # Iteration 1 
        Err1 = []
         
        # Define initial guess
        r_est = 1

        p_T_0 = [r_est,r_est,r_est]
        Initial_guess = p_T_0

        
        ####################################
        #Final TDOA

        TDOA = t_all
        TDOA = np.add(TDOA,Err)
       
        d = np.multiply(c,TDOA)
        ################################### 
        
        i = 1
        
        f = np.zeros((M-1))
        
        dDEC_f = np.zeros((M-1,3))
         
        for i in range(1,M):
            
            f[i-1] = np.linalg.norm( np.array(p_T_0 -  np.array(P)[i][:]   )) - np.linalg.norm(np.array(p_T_0-  np.array(P)[0][:]) )   
        
            dDEC_f[i-1][0] = ((p_T_0[0] -P[i][0])*np.linalg.norm(np.subtract(p_T_0,P[i][:]))**-1 - (p_T_0[0] -P[0][0])*np.linalg.norm(np.subtract(p_T_0,P[0][:]))**-1)
            
            dDEC_f[i-1][1] = (p_T_0[1]-P[i][1])*np.linalg.norm(np.subtract(p_T_0,P[i][:]))**-1 - (p_T_0[1] -P[0][1])*np.linalg.norm(np.subtract(p_T_0,P[0][:]))**-1
            
            dDEC_f[i-1][2] = (p_T_0[2]-P[i][2])*np.linalg.norm(np.subtract(p_T_0,P[i][:]))**-1 - (p_T_0[2] -P[0][2])*np.linalg.norm(np.subtract(p_T_0,P[0][:]))**-1

        x_nonlin = np.dot((np.linalg.pinv(dDEC_f)),np.transpose(d-f)) + p_T_0 
        
        RA = Cart2Sph(x_nonlin[0],x_nonlin[1],x_nonlin[2])[2]
        DEC = Cart2Sph(x_nonlin[0],x_nonlin[1],x_nonlin[2])[1]
        
        if DECs< 90/rad:
            DEC = DEC - 90/rad
        else:
            DEC = 270/rad - DEC 
        
        Fin_Err_All = []
        Fin_Err = abs(Angsepsph(RAs,DECs,RA,DEC))    
        Fin_Err_All.append(Fin_Err)
        
        ####################

        for i in range(Iter):
               
                Initial_guess = x_nonlin
                p_T_0  =  [x_nonlin[0],x_nonlin[1],x_nonlin[2]]
                d = np.multiply(c,TDOA)
                 
                f = np.zeros(M-1)
                dDEC_f = np.zeros((M-1,3))
                 
                for i in range(1,M):
                    
                    f[i-1] = (np.linalg.norm( np.array(p_T_0 -  np.array(P)[i][:]   )) - np.linalg.norm(np.array(p_T_0 -  np.array(P)[0][:]) ))   
                
                    dDEC_f[i-1][0] = ((p_T_0[0] -P[i][0])*np.linalg.norm(np.subtract(p_T_0,P[i][:]))**-1 - (p_T_0[0] -P[0][0])*np.linalg.norm(np.subtract(p_T_0,P[0][:]))**-1)
                    
                    dDEC_f[i-1][1] = (p_T_0[1]-P[i][1])*np.linalg.norm(np.subtract(p_T_0,P[i][:]))**-1 - (p_T_0[1] -P[0][1])*np.linalg.norm(np.subtract(p_T_0,P[0][:]))**-1
                    
                    dDEC_f[i-1][2] = (p_T_0[2]-P[i][2])*np.linalg.norm(np.subtract(p_T_0,P[i][:]))**-1 - (p_T_0[2] -P[0][2])*np.linalg.norm(np.subtract(p_T_0,P[0][:]))**-1
                
                x_nonlin = np.dot((np.linalg.pinv(dDEC_f)),np.transpose(d-f)) + p_T_0 
                
                RA = Cart2Sph(x_nonlin[0],x_nonlin[1],x_nonlin[2])[2]
                DEC = Cart2Sph(x_nonlin[0],x_nonlin[1],x_nonlin[2])[1]
                
            
                if RA<1/rad:
                    RA = 180/rad+RA
 
                if DEC - DECs > 20/rad:
                    DEC = 180/rad - DEC
                
                Fin_Err = abs(Angsepsph(RAs,DECs,RA,DEC))  
                Fin_Err_All.append(Fin_Err)
                
                Err_all.append(Fin_Err)
                RA_all.append(RA)
                DEC_all.append(DEC)
               
        Fin_Err = min(Err_all)
        
        index = np.where(Err_all == min(Err_all))
    
        RA_ind = RA_all[index[0][0]]
        DEC_ind = DEC_all[index[0][0]]
        
        return Fin_Err,RA_ind*rad,DEC_ind*rad,RAs*rad,DECs*rad
 


def Friedlander_Err(P,Err,RA_set,DEC_set):   
    
    Fried_Err_all = []
    
    combination4 = combinations(P,4)
    Com_4  = []
    for i in combination4:
        Com_4.append(i)
     
    for i in range(0,len(Com_4)):
        
        permutation4 =  permutations(Com_4[i],4)
        Perm_4 = []
        for i in permutation4:
            Perm_4.append(i)
        
    for i in range(0,len(Com_4)):
            
            Fried_Err_all.append(Friedlander(Com_4[i],     [Err[0],Err[1],Err[2]],RA_set,DEC_set)[0])
            
    return min(Fried_Err_all)






def Fang_Err(P,Err,RA_set2,DEC_set2):
    Fang_Err_all = []
    
    combination4 = combinations(P,3)
    Com_4  = []
    for i in combination4:
        Com_4.append(i)
       
    for i in range(0,len(Com_4)):
      
        Fang_Err_all.append(FangD3(Com_4[i][0], Com_4[i][1],Com_4[i][2]       ,[Err[0],Err[1]],RA_set2,DEC_set2) [0] )
        Fang_Err_all.append(FangD3(Com_4[i][0], Com_4[i][2],Com_4[i][1]       ,[Err[0],Err[1]],RA_set2,DEC_set2)[0] )
        Fang_Err_all.append(FangD3(Com_4[i][1], Com_4[i][2],Com_4[i][0]       ,[Err[0],Err[1]],RA_set2,DEC_set2)[0] )
        Fang_Err_all.append(FangD3(Com_4[i][1], Com_4[i][0],Com_4[i][2]       ,[Err[0],Err[1]],RA_set2,DEC_set2)[0] )
        Fang_Err_all.append(FangD3(Com_4[i][2], Com_4[i][0],Com_4[i][1]       ,[Err[0],Err[1]],RA_set2,DEC_set2)[0] )
        Fang_Err_all.append(FangD3(Com_4[i][2], Com_4[i][1],Com_4[i][0]       ,[Err[0],Err[1]],RA_set2,DEC_set2)[0] )
        
    Fang_Err_all_2 = []
    
    for i in range(0,len(Fang_Err_all)):
        if type(Fang_Err_all[i]) == np.float64:
            if np.isfinite(Fang_Err_all[i]):
                Fang_Err_all_2.append(Fang_Err_all[i])
        
    return np.abs(min(Fang_Err_all_2 ))




def Err_RA(DEC):

    RA = [DEC]# np.linspace(0,360,1000)#np.linspace(166,180)#np.linspace(0,360,361)
    
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
            
            Err_Foy_all.append(Foy(P_TPSP,Uncert_TPSP,RA[i],DEC[j])[0])
            
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




    
Err_RA1 = Err_RA(15)



#################### TIME THE SCRIPT #############################
toc = time.perf_counter()  
print("Total time = " +str(toc - tic) +"seconds")
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)



























