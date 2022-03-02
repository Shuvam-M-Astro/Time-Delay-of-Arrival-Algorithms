import numpy as np
import matplotlib.pyplot as plt
import math
import sys
np.set_printoptions(threshold=sys.maxsize)
from numpy import loadtxt
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

# Functions to constrain RA and DEC
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

    return Angsep1#,Angsep2

# Function to add more detectors
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

x9 = r_E*np.cos(45/rad) 
z9 = np.sqrt(r_E**2-x9**2)
y9 = 1e-2

x10 = -r_E*np.sin(45/rad)
z10 = np.sqrt(r_E**2-x10**2)
y10 = 1e-3

z11 = -r_E*np.cos(45/rad)
x11 = -1*np.sqrt(r_E**2-z11**2)
y11= 2e-3

z12 = -r_E*np.sin(45/rad)
x12 = np.sqrt(r_E**2-z12**2)
y12 = 3e2

# 156 for equilateral configuration
# 123 for one sided detectors

# 6 Detector Configurations
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

############################################################################
############################################################################
############################################################################
############################################################################


    
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
        ####################  Time-delays from RA and DEC ########################
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

        # Time-delay calculation and add error
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
    

def Fang_Err(P,Err,RA_set2,DEC_set2):
    Fang_Err_all = []
    
    # Combinations of 3 for Fang's algorithm 
    combination4 = combinations(P,3)
    Com_4  = []
    for i in combination4:
        Com_4.append(i)
       
    # Compute for all permutations of the detector configuration and select best one
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


def Uncertainty_Fang(Source,P,Uncert):
    
    Fin = []
    Fin_RA = []
    Fin_DEC = []
    
    if type(Source) == int:
        Source = [Source]
   
    Uncert_Fang = [Uncert,Uncert]
    
    for i in range(0,len(Source)):
      
            RA = Source[0][0]
            DEC = Source[0][1]

            Err_Fang = []

            Uncert_arr = np.arange(0,Uncert,Uncert/1) 

            for m in range(0,len(Uncert_arr)):
                Fang_0 = Fang_Err(P,[Uncert,Uncert],RA,DEC)
            
            Err_Fang.append(Fang_0)
            max_Err =  max(Err_Fang)  
           
   ########################### FANG ###########################################        
            Err = np.arange(0,max_Err,max_Err/50)
   
            for r in range(0,len(Err)):
                for a in range(0,360):
                    Fin_RA.append(RA + Err[r] * np.cos(a))
                    Fin_DEC.append(DEC + Err[r] * np.sin(a))
          
    return max_Err#Fin_RA,Fin_DEC


# Run Fang's algorithm

def Uncertainty_all(Source,P,Uncert):
    
    P1 = np.copy(P)
    Fin_RA = []
    Fin_DEC = []

    RA = Source[0][0]
    DEC = Source[0][1]
    
    Fang = Uncertainty_Fang(Source,P1,Uncert)   
   
    max_Err = Fang
    
    Err = np.arange(0,max_Err,max_Err/50)
   
    RA_min,RA_max = RA-max_Err,RA+max_Err
    DEC_min,DEC_max = DEC-max_Err,DEC+max_Err, 
    
    # Circle with radius as max_Err
    for r in range(0,len(Err)):
        for a in range(0,360):
            Fin_RA.append(RA + Err[r] * np.cos(a))
            Fin_DEC.append(DEC + Err[r] * np.sin(a))

    return Fin_RA,Fin_DEC

    





##############################################################################
##############################################################################
######################         Uncertainty Plot          #####################
##############################################################################
##############################################################################

P = []

# Hexagonal Configuration
P.append(D1)
P.append(D2)
P.append(D3)

P.append(D5)
P.append(D6)
P.append(D7)

Fin_RA = []
Fin_DEC = []

# Selected positions of source for computations
RAar = [165,15,115,65]
DECar = [15,45,75,-15,-45,-75]

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

for i in range(0,len(RAar)):
    for j in range(0,len(DECar)):
        True_S.append([RAar[i],DECar[j]])
for i in range(0,len(True_S)):
    RA_1.append(True_S[i][0])
    DEC_1.append(True_S[i][1])   

ax.scatter(np.divide(RA_1,rad),np.divide(DEC_1,rad),marker='x',color='r',s=0.5,alpha=1,label='True position');  

ax.legend(loc='upper left')
xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
ax.set_xticklabels(xlab, weight=8,fontsize=4)
plt.legend(loc='upper left',fontsize=4,facecolor='white', framealpha=1)
ax.grid(color='k', linestyle='dotted', linewidth=0.5)














