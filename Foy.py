
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

# Set distance to the source
R_set = 1e31

#Define earth radius
r_E = 6.371e6 #m

#Convert radian to degrees
rad = 57.29577951308232

#Define the speed of light
c = 299792458;

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
            
        # Define detector configuration - Equilateral triangle around the equator for now
        
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
    
        R_set = 1e31
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
       
        t12 = np.cos(Ang2s) * d12 / c 
        t12 = t12 + Err[0]
         
        t13 = np.cos(Ang3s) * d13 / c 
        t13 = t13 + Err[1]
         
        t14 = np.cos(Ang4s) * d14 / c 
        t14 = t14 + Err[2]
    
        # Left side and right side of the equation
    
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




def Uncertainty_Friedlander(RA,DEC,P,Width):
    Fin = []
    Fin_RA = []
    Fin_DEC = []
    
    #if len(RA) == 0:
    RA = [RA]
        
    if type(DEC) == int:
        DEC = [DEC]
    
    ########################### FRIED ###########################################         
    
    for i in range(0,len(RA)):
        for j in range(0,len(DEC)):

            Err_Fried = []
            
            Uncert_arr = np.arange(0,Width,Width/10) 

            for m in range(0,len(Uncert_arr)):
                Err_Fried.append(Friedlander_Err(P,[Width,Width,0],RA[i],DEC[j]))

            max_Err =  max(Err_Fried)        
   ########################### FRIED ###########################################         
           
            Err = np.arange(0,max_Err,max_Err/50)

            RA_min,RA_max = RA[i]-max_Err,RA[i]+max_Err
            DEC_min,DEC_max = DEC[j]-max_Err,DEC[j]+max_Err, 

            
            for r in range(0,len(Err)):
                for a in range(0,360):
                    Fin_RA.append(RA[i] + Err[r] * np.cos(a))
                    Fin_DEC.append(DEC[j] + Err[r] * np.sin(a))
           
    return max_Err#Fin_RA,Fin_DEC



def Foy(Pini,Err,RAs,DECs,RAest,DECest):
        
        if DECs<0:
            DECs = 180 + DECs#90 + -1*DECs
    
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

        ###################################################
        ###################################################
        # After time delays
        ###################################################
        ###################################################
        
        # Iteration 1 
        Err1 = []
         
        # Define initial guess
        r_est = 10

        A = Sph2Cart(r_est,DECest/rad,RAest/rad)
        
        p_T_0 = [A[0],A[1],A[2]]
        Initial_guess = p_T_0

        ####################################
        #Add Uncertainty to Time-delay

        TDOA = t_all
     
        Err_0 = np.zeros(len(P)-2)
        for i in range(0,len(Err_0)):
            Err_0[i] = Err
        Err_0 = np.ndarray.tolist(Err_0)
        
        Err_0 = Err_0 + [0]
        
        TDOA = np.add(TDOA,Err_0)
      
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
       
        # In case of 0 Angsep error
        if Err == 0:
            Fin_Err = max(abs(RAs - RA), abs(DECs - DEC))
        else:
            Fin_Err = min(Err_all)
      
        index = np.where(Err_all == min(Err_all))
    
        RA_ind = RA_all[index[0][0]]
        DEC_ind = DEC_all[index[0][0]]
       
        return Fin_Err,RA_ind*rad,DEC_ind*rad,RAs*rad,DECs*rad
 
def Foy_all(Pini,Err,RAs,DECs,RAest,DECest):
    
    combination4 = combinations(Pini,4)
    Com_4  = []
    for i in combination4:
        Com_4.append(i)

    Err_Foy_all = [   ]
    
    for i in range(len(Com_4)):
        Err_Foy_all.append( Foy(Com_4[i]  ,Err,RAs,DECs,RAest,DECest     )[0])
   
    return min((Err_Foy_all))


def Uncertainty_Foy(RA,DEC,P,Width):
    
    Fin = []
    Fin_RA = []
    Fin_DEC = []
    
    RA = [RA]
        
    DEC = [DEC]
    
    for i in range(0,len(RA)):
        for j in range(0,len(DEC)):
    
            Fried1 = Uncertainty_Friedlander(RA[i],DEC[j],P,Width)
          
            Err_Foy_all = []
            Err_Foy_all2 = []
            
            Uncert_arr = [Width] 
            
            for m in range(0,len(Uncert_arr)):
               
                        Err_Foy_all1 = (Foy_all(P,Width,RA[i],DEC[j],RA[i] + Fried1 ,DEC[j] + Fried1))
                        Err_Foy_all2 = (Foy_all(P,Width/10,RA[i],DEC[j],RA[i] + Fried1 ,DEC[j] + Fried1))
                        if Err_Foy_all1>Err_Foy_all2:
                            Err_Foy_all.append(Err_Foy_all1)
                        else:
                            Err_Foy_all.append(360)
 
            max_Err =  min(Err_Foy_all)       
  
            Err = np.arange(0,max_Err,max_Err/50)
   
            RA_min,RA_max = RA[i]-max_Err,RA[i]+max_Err
            DEC_min,DEC_max = DEC[j]-max_Err,DEC[j]+max_Err, 
           
            for r in range(0,len(Err)):
                for a in range(0,360):
                    Fin_RA.append(RA[i] + Err[r] * np.cos(a))
                    Fin_DEC.append(DEC[j] + Err[r] * np.sin(a))
            
    if max_Err > Fried1:
        max_Err = Fried1
    
    
    if max_Err<50:
        return max_Err#Fin_RA,Fin_DEC
    else:
        return 1e2

def Uncertainty_Foy_Fin(RA,DEC,P,Width):
    Uncertainty_Foy_Fin1 = Uncertainty_Foy(RA,DEC,P,Width)
    if Uncertainty_Foy_Fin1<50:
        return Uncertainty_Foy_Fin1
    while Uncertainty_Foy_Fin1>50:
        Width+= 1e-4
        Uncertainty_Foy_Fin1 = Uncertainty_Foy(RA,DEC,P,Width)
        
    return Uncertainty_Foy_Fin1
        
        
        
        
        
        
        
        
        