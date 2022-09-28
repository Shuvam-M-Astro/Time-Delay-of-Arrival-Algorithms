#include<iostream>
#include <stdio.h>     
#include <math.h>   
#include <fstream>
#include <random>
#include <ctime>
#include <cmath>   
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <chrono>

using namespace std;
using namespace std::chrono;

/* Initalize static variables */
#define matsize 3
#define PI 3.14159265
#define rad 57.29577951308232
#define r_E 6.371e6 //m

using namespace std;

const double SMALL = 1.0E-30;          // used to stop divide-by-zero
const double NEARZERO = 1.0e-10;       // helps in printing

using vec    = vector<double>;         // vector
using matrix = vector<vec>;            // matrix

///////////////////////////////////////
///////////////////////////////////////
///// BASIC INDEPENDANT FUNCTIONS ///
///////////////////////////////////////
///////////////////////////////////////

int dot_product(int vector_a[], int vector_b[]) {
   int product = 0;
   for (int i = 0; i < matsize; i++)
   product = product + vector_a[i] * vector_b[i];
   return product;
}

double dot_product2(double vector_a[], double vector_b[]) {
   double product = 0;
   for (int i = 0; i < matsize; i++)
   product = product + vector_a[i] * vector_b[i];
   return product;
}

void cross_product(int vector_a[], int vector_b[], int temp[]) {
   temp[0] = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1];
   temp[1] = -(vector_a[0] * vector_b[2] - vector_a[2] * vector_b[0]);
   temp[2] = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0];
}

double * cross_product2(double vector_a[], double vector_b[], double temp[]) {
   temp[0] = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1];
   temp[1] = -(vector_a[0] * vector_b[2] - vector_a[2] * vector_b[0]);
   temp[2] = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0];

   return temp;
}


void print( const string &title, const matrix &A )
{
   if ( title != "" ) cout << title << '\n';

   for ( auto &row : A )
   {
      for ( auto x : row ) cout << setw( 15 ) << ( abs( x ) < NEARZERO ? 0.0 : x );
      cout << '\n';
   }
}

//======================================================================

matrix matmul( const matrix &A, const matrix &B )          // Matrix times matrix
{
   int rowsA = A.size(),   colsA = A[0].size();
   int rowsB = B.size(),   colsB = B[0].size();
   assert( colsA == rowsB );

   matrix C( rowsA, vec( colsB, 0.0 ) );
   for ( int i = 0; i < rowsA; i++ )
   {
      for ( int j = 0; j < colsB; j++ )
      {
         for ( int k = 0; k < colsA; k++ ) C[i][j] += A[i][k] * B[k][j];
      }
   }
   return C;
}

//======================================================================

matrix subtract( const matrix &A, const matrix &B )        // Subtract matrices
{
   int rows = A.size(),   cols = A[0].size();
   assert( rows == B.size() && cols == B[0].size() );

   matrix result( rows, vec( cols ) );
   for ( int i = 0; i < rows; i++ )
   {
      for ( int j = 0; j < cols; j++ ) result[i][j] = A[i][j] - B[i][j];
   }
   return result;
}

//======================================================================

matrix oppsign( matrix A )                                  // Minus matrix
{
   for ( auto &row : A )
   {
      for ( auto &e : row ) e = -e;
   }
   return A;
}

//======================================================================

matrix subMatrix( const matrix &A, int i1, int i2, int j1, int j2 )
{
   int rows = i2 - i1 + 1, cols = j2 - j1 + 1;
   matrix result( rows, vec( cols ) );
   for ( int i = i1, r = 0; i <= i2; i++, r++ )
   {
      auto it1 = A[i].begin() + j1, it2 = A[i].begin() + j2 + 1;
      copy( it1, it2, result[r].begin() );
   }
   return result;
}

//======================================================================

matrix assembly( const matrix &A11, const matrix &A12, const matrix &A21, const matrix &A22 )
{
   int k = A11.size();           
   int n = k + A22.size();
   matrix result( n, vec( n ) );

   for ( int i = 0; i < k; i++ )
   {
      copy( A11[i].begin(), A11[i].end(), result[i].begin()     );
      copy( A12[i].begin(), A12[i].end(), result[i].begin() + k );
   }

   for ( int i = k; i < n; i++ )
   {
      copy( A21[i-k].begin(), A21[i-k].end(), result[i].begin()     );
      copy( A22[i-k].begin(), A22[i-k].end(), result[i].begin() + k );
   }

   return result;
}

//======================================================================

matrix inverse( const matrix &A )
{
   int n = A.size();
   if ( n == 1 ) 
   {
      double value = A[0][0];
      if ( abs( value ) < SMALL )
      {
         cerr << "Non-invertible. Giving up.\n";
         exit( 0 );
      }
      return matrix( 1, vec( 1, 1.0 / value ) );
   }

   // Partition into four
   int k = n / 2;
   matrix A11 = subMatrix( A, 0, k - 1, 0, k - 1 );
   matrix A12 = subMatrix( A, 0, k - 1, k, n - 1 );
   matrix A21 = subMatrix( A, k, n - 1, 0, k - 1 );
   matrix A22 = subMatrix( A, k, n - 1, k, n - 1 );

   // Strassen steps
   matrix R1  = inverse( A11 );
   matrix R2  = matmul( A21, R1 );
   matrix R3  = matmul( R1, A12 );
   matrix R4  = matmul( A21, R3 );
   matrix R5  = subtract( R4, A22 );
   matrix R6  = inverse( R5 );
   matrix X12 = matmul( R3, R6 );
   matrix X21 = matmul( R6, R2 );
   matrix R7  = matmul( R3, X21 );
   matrix X11 = subtract( R1, R7 );
   matrix X22 = oppsign( R6 );

   return assembly( X11, X12, X21, X22);
}

void multiply(double arr1[][3], double arr2[][3], double arr3[][3]){
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 3; j++){
            arr3[i][j] = 0;
            for(int k = 0; k < 3; k++){
                arr3[i][j] += arr1[i][k] * arr2[k][j];
            }
        }
    }
}

void add(double A[][3], double B[][3], double C[][3])
{
    int i, j;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            C[i][j] = A[i][j] + B[i][j];
}

double norm (double ai, double aj, double ak)
{
  double A = sqrt(pow(ai,2) +pow(aj,2) +  pow(ak,2) );
  return A;
}

void scalarProductMat(double mat[][3], double k)
{
    // scalar element is multiplied by the matrix
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            mat[i][j] = mat[i][j] * k;       
}


std::random_device rd;
std::mt19937 gen(rd());
int random_number(int low, int high)
{
    std::uniform_int_distribution<> dist(low, high);
    return dist(gen);
}
 

int random_choice(int a, int b)
{
      srand(time(NULL));
   
      int r = rand()%2;

      if(r==0)
            return a;
      else
            return b;
}

/* Function to square */

float Square(float value){
    // Multiply value two times
    return value*value;
}

float * Cart2Sph(float x, float y, float z ) {

   float* Vals = new float[5]; 
   if(x == 0 ){
    x = 1e-100;
   }
   float DEC = acos (z / sqrt (x*x + y*y + z*z) ) ;
   float RA = atan(y/x);
   float R = sqrt(x*x + y*y + z*z);
   Vals[0] = R;
   Vals[1] = DEC;
   Vals[2] = RA;

   return Vals;

}

float * Sph2Cart(float R, float El, float Az) {

   float* Vals = new float[5]; 
   float x = R * cos(Az) * sin(El);
   float y = R * sin(Az) * sin(El);
   float z = R * cos(El);
   Vals[0] = x;
   Vals[1] = y;
   Vals[2] = z;

   return Vals;
}

float Angsep(float x1, float y1, float z1, float x2, float y2, float z2, float x3, float y3, float z3){

    float a2  = sqrt( pow((x2-x1),2) + pow((y2-y1),2) + pow((z2-z1),2)              );
    float b2  = sqrt( pow((x3-x1),2) + pow((y3-y1),2) + pow((z3-z1),2)              );
    float c2  = sqrt( pow((x3-x2),2) + pow((y3-y2),2) + pow((z3-z2),2)              );

    float Angsep1 = acos( (pow(a2,2) + pow(b2,2) - pow(c2,2)) / (2 * a2 * b2)  );

    return Angsep1;

}

void transpose(float A[][3], float B[][3])
{
    int i, j;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            B[i][j] = A[j][i];
}

void transpose2(double A[][3], double B[][3])
{
    int i, j;
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++)
            B[i][j] = A[j][i];
}


float Angsepcart(float x1, float y1, float z1, float x2, float y2, float z2){

    float *Vals1;
    Vals1 = Cart2Sph(x1,y1,z1);

    float RA1 = *(Vals1 + 2);
    float Dec1 = PI/2 - *(Vals1 + 1);

    float *Vals2;
    Vals2 = Cart2Sph(x2,y2,z2);

    float RA2 = *(Vals2 + 2);
    float Dec2 = PI/2 - *(Vals2 + 1);

    float Angsep1 = acos( sin(Dec1) * sin(Dec2) + cos(Dec1) * cos(Dec2) * cos(RA1 - RA2) ) * rad;
    return Angsep1;
}
 
float length(float x1, float y1, float z1, float x2, float y2, float z2){
    float A = sqrt( pow((x1-x2),2) + pow((y1-y2),2)     +   pow((z1-z2),2)       );
    return A;
}


// Functions to constrain RA and DEC

float RA(float RA_in){
    if(RA_in > PI*2){
        float RA = RA_in  - 2*PI;
        return RA;
    }
    else if(RA_in < 0){
        float RA = 2*PI + RA_in;
        return RA;
    }
    else if(0 < RA_in < PI*2){
        float RA = RA_in;
        return RA;
    }
    
}

float DEC(float DEC_in){
    if(DEC_in > PI/2){
        float DEC = PI - DEC_in;
        return DEC;
    }
    else if(DEC_in < -PI/2){
        float DEC = -PI - DEC_in;
        return DEC;
    }
    else if(-PI/2 < DEC_in < PI/2){
        float DEC =  DEC_in;
        return DEC;
    }

    
}

int Mat_norm(double Mat[3][3]){
    int r = 3, c = 3;
  
  
    int sum_of_sq = 0;
    for(int i=0; i<r ;i++)
  {
    for(int j=0; j<c; j++)
    {
      sum_of_sq += (Mat[i][j] * Mat[i][j]); 
    }
  }
  
  float out;
  out = sqrt(sum_of_sq);

    return out;
}

float Angsepsph(float RA1, float Dec1,  float RA2, float Dec2){

    double Angsep1 = acos( sin(Dec1)*sin(Dec2) + cos(Dec1)*cos(Dec2)*cos(RA1 - RA2)   )*rad;

    return Angsep1;
}

///////////////////////////////////////
//////////////////////////////////////
///// DEPENDANT FUNCTIONS ///////////
///////////////////////////////////////
//////////////////////////////////////

float more(){

    float N_in = 100;

    float A_all_1000[100] = {};
    float B_all_1000[100] = {};
    float C_all_1000[100] = {};

    for (int input = 1; input < N_in; input++) {
		
        float A = random_number(0, 6.371e6);
        A = random_choice(-A,A);
        float A_left = sqrt( pow(r_E,2)  - pow(A,2)      );
        float B = random_number(0,A_left);
        B = random_choice(B,-B);
        float B_left = sqrt(pow(r_E,2) - pow(A,2) - pow(B,2)     );
        float C = -B_left;
        C = random_choice(-B_left,B_left);

        A_all_1000[input] = A;
        B_all_1000[input] = B;
        C_all_1000[input] = C;

        std::ofstream outfile;
        outfile.open("A.txt", std::ios_base::app); // append instead of overwrite
        outfile << A << "\n" ; 
         
        std::ofstream outfile1;
        outfile1.open("B.txt", std::ios_base::app); // append instead of overwrite
        outfile1 << B << "\n" ; 
       
       std::ofstream outfile2;
        outfile2.open("C.txt", std::ios_base::app); // append instead of overwrite
        outfile2 << C << "\n" ; 
    }

    return 0;
}


float more_eq(){

    float N_in = 100;

    float A_all_1000[100] = {};
    float B_all_1000[100] = {};
    float C_all_1000[100] = {};

    for (int input = 1; input < N_in; input++) {
		
        float A = random_number(0, 6.371e6);
        A = random_choice(-A,A);
        float A_left = sqrt( pow(r_E,2)  - pow(A,2)      );
        float B = random_number(0,A_left);
        B = random_choice(B,-B);
        float B_left = sqrt(pow(r_E,2) - pow(A,2) - pow(B,2)     );
        float C = 1e-100;
    

        A_all_1000[input] = A;
        B_all_1000[input] = B;
        C_all_1000[input] = C;

        std::ofstream outfile;
        outfile.open("A.txt", std::ios_base::app); // append instead of overwrite
        outfile << A << "\n" ; 
         
        std::ofstream outfile1;
        outfile1.open("B.txt", std::ios_base::app); // append instead of overwrite
        outfile1 << B << "\n" ; 
       
       std::ofstream outfile2;
        outfile2.open("C.txt", std::ios_base::app); // append instead of overwrite
        outfile2 << C << "\n" ; 
    }

    return 0;
}

 
void * rot(double ai, double aj, double ak, double bi, double bj, double bk, double Fin_Add[][3]){

    double a[3] = {ai,aj,ak};

    double a_fin[] = {ai/ norm(ai,aj,ak), aj/ norm(ai,aj,ak)    ,  ak/ norm(ai,aj,ak)    };
   
    double b[3] = {bi,bj,bk};
    double b_fin[] = {bi/ norm(bi,bj,bk), bj/ norm(bi,bj,bk)    ,  bk/ norm(bi,bj,bk)    };
    
    double A[] = {a_fin[0],a_fin[1],a_fin[2]};
    double B[] = {b_fin[0],b_fin[1],b_fin[2]};
   
    double temp[3];

    double * Cross;

    Cross = cross_product2(A,B,temp);
    
    double Dot = dot_product2(A,B);

    double s = norm(Cross[0],Cross[1],Cross[2]);
    
    double kmat[3][3] = {{0,-Cross[2],Cross[1]},{Cross[2],0,-Cross[0]},{Cross[1],Cross[0],0}};
   
    double Eye[3][3] = {{1,0,0},{0,1,0},{0,0,1}};

    double Add[3][3];

    add(Eye,kmat,Add);

    double kmat2[3][3];
    multiply(kmat,kmat,kmat2);

    double Fin_mul = (1 - Dot) / (pow(s,2));

    scalarProductMat(kmat2,Fin_mul);
    
    double Fin_add[3][3];

    add(Add,kmat2, Fin_add); 

    return 0;
}

float isRotationMatrix(double R[3][3]){

    double Rt[3][3], i, j;
    transpose2(R, Rt);
    double shouldbeIdentity[3][3];
    multiply(Rt, R,shouldbeIdentity);
    
    // Create Identity matrix
    double mat2[3][3] = {
        1,0,0,
        0,1,0,
        0,0,1,
    };
    
    double mat_n[3][3];

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            mat_n[i][j] = (mat2[i][j] - shouldbeIdentity[i][j]);	
        }
    }

    double norm_mat_n = Mat_norm(mat_n) ;
    
    return norm_mat_n;

}


double rotationMatrixToEulerAngles(double R[3][3], double Ang[3]){


    double B = isRotationMatrix(R);

    double sy = sqrt( R[0][0] * R[0][0] + R[1][0] * R[1][0]      );
    
    if (sy > 1e-6){
        double x = atan2(R[2][1],R[2][2]);
        double y = atan2(-R[2][0], sy);
        double z = atan2(R[1][0], R[0][0]);

        Ang[0] = x;
        Ang[1] = y;
        Ang[2] = z;
        
    }
    else if (sy < 1e-6){
      
        double x = atan2(-R[1][2], R[1][1]);
        double y = atan2(-R[2][0], sy);
        double z = 0;
        
        Ang[0] = x;
        Ang[1] = y;
        Ang[2] = z;
        
    }

    return 0;
    
    
}



void combinationUtil(float arr[], float data[],
                    int start, int end,
                    int index, int r)
{

    if (index == r)
    {
        for (int j = 0; j < r; j++)
            cout << data[j] << " ";
        cout << endl;
        return;
    }
 
    for (int i = start; i <= end &&
        end - i + 1 >= r - index; i++)
    {
        data[index] = arr[i];
        combinationUtil(arr, data, i+1,
                        end, index+1, r);
    }
}



/***************************
 * **********************
 * *************************/

float FangD3(float D1[3], float D2[3], float D3[3], float Err, float RA_set, float DEC_set){
    /* DEC_set constrainst  */
    if (DEC_set<90){
        DEC_set == 90 - DEC_set;}

    /*
     x,y,z coordinates of the detector config  */
    float x1 = D1[0];
    float y1 = D1[0];
    float z1 = D1[0];

    float x2 = D2[0];
    float y2 = D2[0];
    float z2 = D2[0];

    float x3 = D3[0];
    float y3 = D3[0];
    float z3 = D3[0];

    float x1_1 = x1;
    float y1_1 = y1;
    float z1_1 = z1;

    /* Translate detector configuration*/

    x1 = x1 - x1_1 + 1e-6;
    y1 = y1 - y1_1;
    z1 = z1 - z1_1;
   
    x2 = x2 - x1_1 ;
    y2 = y2 - y1_1;
    z2 = z2 - z1_1;

    x3 = x3 - x1_1 ;
    y3 = y3 - y1_1;
    z3 = z3 - z1_1;

    float X1[3] = {x1,y1,z1};
    float X2[3] = {x2,y2,z2};
    float X3[3] = {x3,y3,z3};

     /* Source coorinates  */

    float P[3][3] = {{X1[1],X1[2],X1[3]},   {X2[1],X2[2],X2[3]},    {X3[1],X3[2],X3[3]}};

    float DECs = DEC_set;
    float RAs = RA_set/rad;

    /*calculate spherical coordinates of source from cartesian
    same with detectors */

    float DEC_all[3] = {0,0,0};
    float RA_all[3] = {0,0,0};  

     for (int i : {0, 1, 2}){
        DEC_all[i] == DEC(PI/2 - Cart2Sph(P[i][0],P[i][1],P[i][2])[1] );
        RA_all[i] == RA(Cart2Sph(P[i][0],P[i][1],P[i][2])[2]);

     ;}

    float DECar = PI/2 - DECs/rad;
    float RAar  = RAs/rad;
    DECar == DEC(DECar);
    RAar == RA(RAs);

    float DEC1 = PI/2 - Cart2Sph(x1,y1,z1)[1]/rad;
    float RA1 = Cart2Sph(x1,y1,z1)[2] / rad;
    DEC1 == DEC(DEC1);
    RA1 == RA(RA1);

    float DEC2 = PI/2 - Cart2Sph(x2,y2,z2)[1]/rad;
    float RA2 = Cart2Sph(x2,y2,z2)[2] / rad;
    DEC2 == DEC(DEC2);
    RA2 == RA(RA2);

    float DEC3 = PI/2 - Cart2Sph(x3,y3,z3)[1]/rad;
    float RA3 = Cart2Sph(x3,y3,z3)[2] / rad;
    DEC3 == DEC(DEC3);
    RA3 == RA(RA3);

    /*   Spherical coordinates of origin   */
    float DEC0 = 0;
    float RA0 = 0;

    /*  Define distances between the detectors   */

    float d21 = sqrt( pow((x2-x1),2) + pow((y2-y1),2) + pow((z2-z1),2)                     );
    float d31 = sqrt( pow((x3-x1),2) + pow((y3-y1),2) + pow((z3-z1),2)                     );
    float d12 = d21;
    float d13 = d31;

    float r_all[2][3] = {{x2-x1,y2-y1,z2-z1},{x3-x1,y3-y1,z3-z1}};

    float r11 = r_all[0][0];
    float r12 = r_all[0][1];
    float r13 = r_all[0][2];

    float r21 = r_all[1][0];
    float r22 = r_all[1][1];
    float r23 = r_all[1][2];

    float ri = r12*r23 - r22*r13;
    float rj = r11*r23 - r21*r13;
    float rk = r11*r22 - r21*r12;

    float a = ri;
    float b = rj;
    float c = rk;
    float d = ri*x1 + y1*rj + z1*rk;

    /* Find the rotation matrix using the normal vectors of the planes  */

    double vector_1[3] = {a,b,c};
    double vector_2[3] = {0,1,0};
    double unit_vector_1[3] = {{vector_1[0] / norm(vector_1[0],vector_1[1],vector_1[2])}, {vector_1[1] / norm(vector_1[0],vector_1[1],vector_1[2])}   , {vector_1[2] / norm(vector_1[0],vector_1[1],vector_1[2])}   };

    double unit_vector_2[3] = {{vector_2[0] / norm(vector_2[0],vector_2[1],vector_2[2])}, {vector_2[1] / norm(vector_2[0],vector_2[1],vector_2[2])}   , {vector_2[2] / norm(vector_2[0],vector_2[1],vector_2[2])}   };

    double dot_product = dot_product2(vector_1,vector_2);
    double angle = acos(dot_product);

    double temp[3] = {};
    cross_product2(vector_1,vector_2,temp);

    double axis[3] = {{temp[0]/norm(temp[0],temp[1],temp[2])}, {temp[1]/norm(temp[0],temp[1],temp[2])}  ,{temp[2]/norm(temp[0],temp[1],temp[2])}};
 
    double cos_angle = cos(angle);
    double s = sqrt(1 - cos_angle*cos_angle);
    double C = 1 - cos_angle;

    double x = axis[0];
    double y = axis[1];
    double z = axis[2];

    /*  Matrix elements for the rotation matrix  */
    double r_mat_11 = x*x*C + cos_angle;
    double r_mat_12 = x*y*C - z*s;
    double r_mat_13 = x*z*C + y*s;

    double r_mat_21 = y*x*C + z*s;
    double r_mat_22 = y*y*C + cos_angle;
    double r_mat_23 = y*z*C - x*s;

    double r_mat_31 = z*x*C - y*s;
    double r_mat_32 = z*y*C + x*s;
    double r_mat_33 = z*z*C + cos_angle;

    /* Final rotation matrix  */
    matrix rmat = {{r_mat_11,r_mat_12,r_mat_13},{r_mat_21,r_mat_22,r_mat_23},{r_mat_31,r_mat_32,r_mat_33}};

    /* Detector coordinates in new frame  */
    /* Rotation */

    matrix rmatinv = inverse(rmat);

    matrix X1mat  = {{X1[0]},{X1[1]},{X1[2]}};
    matrix X1_2x = matmul(rmatinv,X1mat);

    matrix X2mat  = {{X2[0]},{X2[1]},{X2[2]}};
    matrix X2_2x = matmul(rmatinv,X2mat);

    matrix X3mat  = {{X3[0]},{X3[1]},{X3[2]}};
    matrix X3_2x = matmul(rmatinv,X3mat);
 
    std::vector<double> input({ X1_2x[0] });
    float X1_2x2[1];
    std::copy(input.begin(), input.end(), X1_2x2);

    std::vector<double> input1({ X1_2x[1] });
    float X1_2y2[1];
    std::copy(input1.begin(), input1.end(), X1_2y2);

    std::vector<double> input2({ X1_2x[2] });
    float X1_2z2[1];
    std::copy(input2.begin(), input2.end(), X1_2z2);


    std::vector<double> input3({ X2_2x[0] });
    float X2_2x2[1];
    std::copy(input3.begin(), input3.end(), X2_2x2);

    std::vector<double> input4({ X2_2x[1] });
    float X2_2y2[1];
    std::copy(input4.begin(), input4.end(), X2_2y2);

    std::vector<double> input5({ X2_2x[2] });
    float X2_2z2[1];
    std::copy(input5.begin(), input5.end(), X2_2z2);

    std::vector<double> input6({ X3_2x[0] });
    float X3_2x2[1];
    std::copy(input6.begin(), input6.end(), X3_2x2);

    std::vector<double> input7({ X3_2x[1] });
    float X3_2y2[1];
    std::copy(input7.begin(), input7.end(), X3_2y2);

    std::vector<double> input8({ X3_2x[2] });
    float X3_2z2[1];
    std::copy(input8.begin(), input8.end(), X3_2z2);

    /*    Convert to spherical coordinates */
    float DEC1_2 = PI/2 - Cart2Sph(X1_2x2[0],X1_2y2[0],X1_2z2[0])[1];
    float DEC2_2 = PI/2 - Cart2Sph(X2_2x2[0],X2_2y2[0],X2_2z2[0])[1];
    float DEC3_2 = PI/2 - Cart2Sph(X3_2x2[0],X3_2y2[0],X3_2z2[0])[1];

    float RA1_2 = Cart2Sph(X1_2x2[0],X1_2y2[0],X1_2z2[0])[2];
    float RA2_2 = Cart2Sph(X2_2x2[0],X2_2y2[0],X2_2z2[0])[2];
    float RA3_2 = Cart2Sph(X3_2x2[0],X3_2y2[0],X3_2z2[0])[2];

    DEC1_2 == DEC(DEC1_2);
    DEC2_2 == DEC(DEC2_2);
    DEC3_2 == DEC(DEC3_2);

    RA1_2 = RA(RA1_2);
    RA2_2 = RA(RA2_2);
    RA3_2 = RA(RA3_2);

    /*    Compare angles in the new frame to those in   the original frame */

    float d12_2 = length(X1_2x2[0],X1_2y2[0],X1_2z2[0],X2_2x2[0],X2_2y2[0],X2_2z2[0]        );
    float d13_2 = length(X1_2x2[0],X1_2y2[0],X1_2z2[0],X3_2x2[0],X3_2y2[0],X3_2z2[0]        );

   /*Rotation matrices to direct euler angles */
    /* Convert rotation matrix to floats */

    std::vector<double> rmat_11({ rmat[0][0] });
    double rmat11[1];
    std::copy(rmat_11.begin(), rmat_11.end(), rmat11);

    std::vector<double> rmat_12({ rmat[0][1] });
    double rmat12[1];
    std::copy(rmat_12.begin(), rmat_12.end(), rmat12);

    std::vector<double> rmat_13({ rmat[0][2] });
    double rmat13[1];
    std::copy(rmat_13.begin(), rmat_13.end(), rmat13);

    std::vector<double> rmat_21({ rmat[1][0] });
    double rmat21[1];
    std::copy(rmat_21.begin(), rmat_21.end(), rmat21);

    std::vector<double> rmat_22({ rmat[1][1] });
    double rmat22[1];
    std::copy(rmat_22.begin(), rmat_22.end(), rmat22);

    std::vector<double> rmat_23({ rmat[1][2] });
    double rmat23[1];
    std::copy(rmat_23.begin(), rmat_23.end(), rmat23);

    std::vector<double> rmat_31({ rmat[0][2] });
    double rmat31[1];
    std::copy(rmat_31.begin(), rmat_31.end(), rmat31);

    std::vector<double> rmat_32({ rmat[1][2] });
    double rmat32[1];
    std::copy(rmat_32.begin(), rmat_32.end(), rmat32);

    std::vector<double> rmat_33({ rmat[2][2] });
    double rmat33[1];
    std::copy(rmat_33.begin(), rmat_33.end(), rmat33);

    double rmat2[3][3] = {{rmat11[0],rmat12[0],rmat13[0]},{rmat21[0],rmat22[0],rmat23[0]},{rmat31[0],rmat32[0],rmat33[0]}};

    double Ang[3] = {};

    float E_1 = rotationMatrixToEulerAngles(rmat2,Ang);
    float RAar_2 = RAar + Ang[1];
    float DECar_2 = DECar + Ang[2];

   /*    Time-delays from RA  and DEC*/
    float R_set = 1e31;
    float xs = Sph2Cart(R_set,DECar_2,RAar_2)[0];
    float ys = Sph2Cart(R_set,DECar_2,RAar_2)[1];
    float zs = Sph2Cart(R_set,DECar_2,RAar_2)[2];

    float xar = xs;
    float yar = ys;
    float zar = zs;

    float R_set_2 = 1e6;

    float x1_2 = X1_2x2[0];
    float y1_2 = X1_2y2[0];
    float z1_2 = X1_2z2[0];

    float x2_2 = X2_2x2[0];
    float y2_2 = X2_2y2[0];
    float z2_2 = X2_2z2[0];

    float x3_2 = X3_2x2[0];
    float y3_2 = X3_2y2[0];
    float z3_2 = X3_2z2[0];

    float d1ar = sqrt( pow((xar - x1_2),2) + pow((yar - y1_2),2) +  pow((zar - z1_2),2)           );
    float d2ar = sqrt( pow((xar - x2_2),2) + pow((yar - y2_2),2) +  pow((zar - z2_2),2)           );
    float d3ar = sqrt( pow((xar - x3_2),2) + pow((yar - y3_2),2) +  pow((zar - z3_2),2)           );

    float c_2 = 299792458;

    /*    Time Delay calculation and add error*/
    float t21ar = (d2ar - d1ar)/c;
    float t21 = t21ar;
    float t12 = -t21 + Err;

    float t31ar = (d3ar - d1ar)/c;
    float t31 = t31ar;
    float t13 = -t31 + Err;

    /*    Angles between pairs of detectors w.r.t D1*/
    float Br12 = acos(299792458*t12/d12)*rad ;
    float Br13 = acos(299792458*t13/d13)*rad ;
    /*    Find the cource coordinates in new frame*/
    float DDECs_2 = asin(( (cos(Br12/rad) - (cos(DEC2_2) / cos(DEC3_2)) * cos(Br13/rad)) / (sin(DEC2_2) - (cos(DEC2_2) / cos(DEC3_2)) * sin(DEC3_2))                                   ));

    float RAs_2 = RA2_2 + acos(((cos(Br12/rad) - sin(DEC2_2) * sin(DDECs_2)) / (cos(DEC2_2) * cos(DDECs_2))));

    float RAs_2_2 = PI*2 - RAs_2;

    /*    Errors in 2nd frame*/
    float RA2_1_er = (abs((RAs_2*rad  - RAar_2*rad)/(RAar_2*rad))*100);
    float RA2_2_er = (abs((PI*2*rad - RAs_2*rad  - RAar_2*rad)/(RAar_2*rad))*100);
    float DEC2_er = (abs((DDECs_2*rad - DECar_2*rad)/(DECar_2*rad))*100);

    /*    Rotation to original frame of reference*/
    float FinRA_est =  180 -  RAs_2*rad - Ang[1]*rad  ;
    float FinDEC_est = PI/2*rad -( 1 * (DDECs_2*rad - Ang[2]*rad)) ;

    /*    Constrain calculated esetimate*/
    if (FinDEC_est < PI/2*rad){  
            FinDEC_est = FinDEC_est;}
    if (FinDEC_est > PI/2*rad){
        FinDEC_est = PI/2*rad -( 1 * (DDECs_2*rad - Ang[2]*rad)) - PI*rad;}

    if (DECar_2<0){
            DECar_2 = -1*DECar_2;}
        
    if (RAar_2>180/rad){
        FinRA_est = 360 - FinRA_est;}
        
    if (DEC_set<90){
        FinRA_est = 180 - FinRA_est;}
    
    float Fin_Err = abs(Angsepsph(FinRA_est/rad,RAar_2,FinDEC_est/rad,DECar_2));
    
    if (DECar_2*rad<90){
        FinDEC_est = 90 - FinDEC_est;}
        
    if (DECar_2>90/rad){
        FinDEC_est = 0;}

    if (FinRA_est > 360){
        FinRA_est  = FinRA_est - 360;}
    
    if (FinRA_est < 0){
        FinRA_est  = FinRA_est + 360*2;}

    /* In Case the rotation matrix was uninvertible   */

    return Fin_Err;
}


float P_TPSP_LVLH_sDEC_0[3][6] = { {0.127*r_E , 0.048*r_E , -0.041*r_E ,  0.030*r_E, -0.052*r_E, -0.136*r_E}, 
{0.272*r_E ,0.295*r_E, 0.299*r_E ,-0.296*r_E, -0.292*r_E, -0.267*r_E},{
 1e-32,  1e-3,  1e-21 , 1e-12 , 1e-34 , 1e-54} };


float P_TPSP_LVLH_opt_0[3][6] = { {0.127*r_E , 0.048*r_E , -0.041*r_E  , 0.184*r_E, -0.132*r_E, -0.211*r_E},{
 0.272*r_E, 0.295*r_E, 0.299*r_E ,-0.237*r_E ,-0.269*r_E, -0.214*r_E},{
 1e-24 , 1e-6 , 1e-4 , 1e-23  ,1e-34 , 1e-54} };

int main()
{
    
    float D1[3] = {3,4,5};
    float D2[3] = {33,432,5524};
    float D3[3] = {34,45,515};
    float Err = 4;
    float RA_set = 15;
    float DEC_set = 15;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    float DX = FangD3(D1, D2, D3, 1e-3, 12, 12);
    
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count() << "[ns]" << std::endl;
    cout<<12;
    
    return 0;
    
}










