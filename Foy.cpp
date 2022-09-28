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

void subtract2( double A[], double B[], double result[] )        // Subtract matrices
{
   double N = 3;
  
   for ( int i = 0; i < N; i++ )
        result[i] = A[i] - B[i];
   
}

void add2( double A[], double B[], double result[] )        // Subtract matrices
{
   double N = 3;
  
   for ( int i = 0; i < N; i++ )
        result[i] = A[i] + B[i];
   
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

float  Sph2Cart(float R, float El, float Az, float Vals[3]) {

   float x = R * cos(Az) * sin(El);
   float y = R * sin(Az) * sin(El);
   float z = R * cos(El);
   Vals[0] = x;
   Vals[1] = y;
   Vals[2] = z;

   return 0;
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

double Mat_norm(double Mat[3][3]){
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
		//std::cout << "\n" << input;
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
        outfile.open("A.txt", std::ios_base::app); 
        outfile << A << "\n" ; 
         
        std::ofstream outfile1;
        outfile1.open("B.txt", std::ios_base::app); 
        outfile1 << B << "\n" ; 
       
       std::ofstream outfile2;
        outfile2.open("C.txt", std::ios_base::app); 
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
		//std::cout << "\n" << input;
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
        outfile.open("A.txt", std::ios_base::app); 
        outfile << A << "\n" ; 
         
        std::ofstream outfile1;
        outfile1.open("B.txt", std::ios_base::app); 
        outfile1 << B << "\n" ; 
       
       std::ofstream outfile2;
        outfile2.open("C.txt", std::ios_base::app); 
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


void lin_solve(double storedA[3][3], double storedb[3], int n, double x[3]){


//This grand if statement finds and prints out the condition numbers at one and infinity for matrices A with dimensions 2x2.
if(n == 2){
	//Variable cond1 is the condition number at one, and condinf is the condition number at infinity.
	double cond1, condinf;
	//detA2x2 is the determinant of 2x2 matrix, in order to be able to find the inverse. 
	double detA2x2 = storedA[0][0] * storedA[1][1] - storedA[0][1] * storedA[1][0];
	
	//Creates the dynamically allocated inverse matrix invA2x2, since it will be used while computing the condition numbers.
	double** invA2x2 = new double*[2];
	for(int i = 0; i < 2; i++)
	invA2x2[i] = new double[2];
	//Sets the elements of inverse matrix same as the original, only for initialization.
	for(int i = 0; i < 2; i++){
		for(int j = 0; j < 2; j++){
			invA2x2[i][j] = storedA[i][j];
		}
	}
	
	//Finds the inverse of matrix A. 
	double temp;
	temp = invA2x2[0][0];
	invA2x2[0][0] = invA2x2[1][1] / detA2x2;
	invA2x2[0][1] = -invA2x2[0][1] / detA2x2;
	invA2x2[1][0] = -invA2x2[1][0] / detA2x2;
	invA2x2[1][1] = temp / detA2x2;

	//The variables under this row and the for loop finds the condition number at one.
	double maxcolA, maxcolinvA, sumcolA, sumcolinvA = 0;
	for(int i = 0; i < 2; ++i){
		sumcolA = 0;
		sumcolinvA = 0;
		//This for loop finds column sums for both matrix A and its inverse.
		for (int j = 0; j < 2; ++j){
  		    	sumcolA = sumcolA +  fabs(storedA[j][i]);
  		   	sumcolinvA = sumcolinvA + fabs(invA2x2[j][i]);
		}
		//This if statement finds the maximum column sum of A.
		if (sumcolA > maxcolA){
			maxcolA = sumcolA;
		}
		//This if statement finds the maximum column sum of inverse of A.
    		if (sumcolinvA > maxcolinvA){
		     	maxcolinvA = sumcolinvA;
  		}	
	}
	//Condition number at one is the product of maximum column sum of A and its inverse's.
	cond1 = maxcolA * maxcolinvA;
	//Prints out the condition number at one.
	cout << "Condition number at one is: " << cond1 << "\n";

	//The variables under this row and the for loop finds the condition number at infinity.
	double maxrowA, maxrowinvA, sumrowA, sumrowinvA = 0;
	for(int i = 0; i < 2; ++i){
		sumrowA = 0;
		sumrowinvA = 0;
		//This for loop finds row sums for both matrix A and its inverse.
		for (int j = 0; j < 2; ++j){
  		    	sumrowA = sumrowA +  fabs(storedA[i][j]);
  		    	sumrowinvA = sumrowinvA + fabs(invA2x2[i][j]);
		}
		//This if statement finds the maximum row sum of A.
		if (sumrowA > maxrowA){
			maxrowA = sumrowA;
		}
		//This if statement finds the maximum row sum of inverse of A.		
  	  	if (sumrowinvA > maxrowinvA){
 		     	maxrowinvA = sumrowinvA;
  	  	}
	}
	//Condition number at infinity is the product of maximum row sum of A and its inverse's.
	condinf = maxrowA * maxrowinvA;
	//Prints out the condition number at infinity.
	cout << "Condition number at infinity is: " << condinf << "\n" ;
}
//Partial pivoting begins here.
double compA, tempA, tempb, factorA;
//This grand for loop turns matrix A into upper triangular form, using partial pivoting.
for(int i = 0; i < n; i++){
	//Initializes maximum value of pivot, compA and its index, q.
      compA = storedA[i][i];
      int q = i;
      	//This for loop compares each pivot candidate and finds the maximum of them, which becomes compA.
            for(int j = i + 1; j < n; j++)
      	      if(fabs(compA) < fabs(storedA[j][i])){
      			compA = storedA[j][i] ;
            		q = j;
           		}
	//This if statement checks if A is singular and if so, quits.
	//Matrix A is singular if a principal diagonalelement is equal to zero.
	//For machine precision, 10^-(-5) is picked. The values smaller than this number are considered to be zero.
	//if(fabs(compA) < 0.00001){ 
	//	cout << "A is singular.";
	//	return 3;
	//}
	//Swaps the maximum row with the current row for A.
	for(int k = 0; k < n; k++)
      {
            tempA = storedA[q][k];
            storedA[q][k] = storedA[i][k];
            storedA[i][k] = tempA;       
      }
      //Swaps the same rows of b.
      tempb = storedb[q];
      storedb[q] = storedb[i];
      storedb[i] = tempb;
	//Turns A into an upper triangular matrix.
      for(int l = i+1; l < n; l++){
            //factorA sets the l'th element of A equal to zero.
		factorA = storedA[l][i]/storedA[i][i];
            storedb[l] = storedb[l] - factorA * storedb[i];
            //This part sets the elements under the diagonal equal to zero.
            for(int m = 0; m < n; m++){
            	storedA[l][m] = storedA[l][m] - factorA * storedA[i][m];
 			if(fabs(storedA[l][m]) < 0.00001)
 				storedA[l][m] = 0;
 		}
      }
}
//Creates solution vector x with dimensions nx1.

//Starts from bottom for backward substitution.
for(int i = n-1; i >= 0; i--){
	//Sets i'th element of x is equal to the i'th element of b divided by the i'th element of A. 
	x[i] = storedb[i] / storedA[i][i];
	//Substracts the product of x[i] and storedA[j][i] from storedb[j] and finds the solution vector x.
	for(int j = 0; j < i; j++){
		storedb[j] = storedb[j] - storedA[j][i]*x[i];
	}
}

//Prints out the solution vector x.
cout << "The solution vector x is: \n";
for (int i = 0; i < n; i++){
	cout << x[i] << "\n";
}


}


void combinationUtil(float arr[], float data[],
                    int start, int end,
                    int index, int r)
{
    // Current combination is ready
    // to be printed, print it
    if (index == r)
    {
        for (int j = 0; j < r; j++)
            cout << data[j] << " ";
        cout << endl;
        return;
    }
 
    // replace index with all possible
    // elements. The condition "end-i+1 >= r-index"
    // makes sure that including one element
    // at index will make a combination with
    // remaining elements at remaining positions
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

void Friedlander(float P[4],float Err,double RA_set,float DEC_set,float result) {

float RA_set2 = RA_set/rad;
float DEC_set2 = DEC_set/rad;
    
float D1[3] = {P[0]};
float D2[3] = {P[1]};
float D3[3] = {P[2]};
float D4[3] = {P[3]};

float x1 = D1[0];
float y1 = D1[1];
float z1 = D1[2];

float x2 = D2[0];
float y2 = D2[1];
float z2 = D2[2];

float x3 = D3[0];
float y3 = D3[1];
float z3 = D3[2];

float x4 = D4[0];
float y4 = D4[1];
float z4 = D4[2];

float x1_1 = x1;
float y1_1 = y1;
float z1_1 = z1;

/*# Translate detector configuration*/
        
x1 = x1 - x1_1 + 1e-100 ;
y1 = y1 - y1_1  ;
z1 = z1 - z1_1;

x2 = x2 - x1_1 ;
y2 = y2 - y1_1;
z2 = z2 - z1_1;

x3 = x3 - x1_1;
y3 = y3 - y1_1;
z3 = z3 - z1_1;

x4 = x4 - x1_1;
y4 = y4 - y1_1;
z4 = z4 - z1_1;

float R_set = 1e21;
float Vals[3] = {0,0,0};
float Source = {Sph2Cart(R_set, DEC_set, RA_set, Vals)};
float xs = Vals[0];
float ys = Vals[1];
float zs = Vals[2];

xs = xs - x1_1;
ys = ys - y1_1;
zs = zs - z1_1;

/*#Define the distances between the detectors*/
float d21 = sqrt(pow((x2-x1),2)+pow((y2-y1),2)+pow((z2-z1),2));
float d31 = sqrt(pow((x3-x1),2)+pow((y3-y1),2)+pow((z3-z1),2));
float d41 = sqrt(pow((x4-x1),2)+pow((y4-y1),2)+pow((z4-z1),2));

float d12 = d21;
float d13 = d31;
float d14 = d41;

/*#Calculate Spherical coordinates from cartesian coordinates ... Only for angular distance
*/
float RAs = RA_set;
float DECs = DEC_set;

float RAar = RAs;
float DECar = PI/2 - DECs;
DECar = DEC(DECar);
RAar = RA(RAar);

float R1  = Cart2Sph(x1,y1,z1)[0];
float DEC1 = PI/2 - Cart2Sph(x1,y1,z1)[1];
float RA1 = Cart2Sph(x1,y1,z1)[2];
DEC1 = DEC(DEC1);
RA1 = RA(RA1);

float R2  = Cart2Sph(x2,y2,z2)[0];
float DEC2 = PI/2 - Cart2Sph(x2,y2,z2)[1];
float RA2 = Cart2Sph(x2,y2,z2)[2];
DEC2 = DEC(DEC2) ;
RA2 = RA(RA2);

float R3  = Cart2Sph(x3,y3,z3)[0];
float DEC3 = PI/2 - Cart2Sph(x3,y3,z3)[1];
float RA3 = Cart2Sph(x3,y3,z3)[2];
DEC3 = DEC(DEC3);
RA3 = RA(RA3);

float R4  =           Cart2Sph(x4,y4,z4)[0];
float DEC4 = PI/2 - Cart2Sph(x4,y4,z4)[1];
float RA4 =           Cart2Sph(x4,y4,z4)[2];
DEC4 = DEC(DEC4);
RA4 = RA(RA4);

/*Calculate time-dDECays from RA and DEC*/

float xar = xs;
float yar = ys;
float zar = zs;

float d1ar = sqrt( pow((xar - x1),2) + pow((yar - y1),2) + pow((zar - z1),2)   );
float d2ar = sqrt( pow((xar - x2),2) + pow((yar - y2),2) + pow((zar - z2),2)   );
float d3ar = sqrt( pow((xar - x3),2) + pow((yar - y3),2) + pow((zar - z3),2)   );
float d4ar = sqrt( pow((xar - x4),2) + pow((yar - y4),2) + pow((zar - z4),2)   );

float c = 299792458 ;

float t21ar = (d2ar-d1ar)/c ;
float  t21 = t21ar; 
float  t12 = -t21 + Err;

float t31ar = (d3ar-d1ar)/c ;
float  t31 = t31ar; 
float  t13 = -t31 + Err;

float t41ar = (d4ar-d1ar)/c ;
float  t41 = t41ar; 
float  t14 = -t41 + Err;

float Ang2s = Angsepsph(RA2,DEC2,RAar,DECar)/rad;
float Ang3s = Angsepsph(RA3,DEC3,RAar,DECar)/rad;
float Ang4s = Angsepsph(RA4,DEC4,RAar,DECar)/rad;

/* Time-delay and add uncertainty  */
t12 == cos(Ang2s) * d12 / c ;
   
t12 = t12 + Err;
    
t13 == cos(Ang3s) * d13 / c ;
t13 = t13 + Err;
    
t14 == cos(Ang4s) * d14 / c ;
t14 = t14 + Err;

/* Left side and right side  */

float S_12_Z[3] = {{x1-x2},{y1-y2},{z1-z2}};
float S_13_Z[3] = {{x1-x3},{y1-y3},{z1-z3}};
float S_14_Z[3] = {{x1-x4},{y1-y4},{z1-z4}};

double L[3][3] = {{S_12_Z[0],S_12_Z[1],S_12_Z[2]},{S_13_Z[0],S_13_Z[1],S_13_Z[2]},{S_14_Z[0],S_14_Z[1],S_14_Z[2]}}; 

double R[3] = { d12*t12/(d12/c),d13*t13/(d13/c),d14*t14/(d14/c)    };

int n = 3;
double U[3] = {};

lin_solve( L, R, n, U);

float RA_calc = Cart2Sph(U[0],U[1],U[2])[2]*rad;

if (RA_calc<0){
            RA_calc = RA_calc + PI*rad;}

float DEC_calc = Cart2Sph(U[0],U[1],U[2])[1]*rad;

if (RAar > 180/rad){
            RA_calc = 180+RA_calc;}
        
if (DEC_set>0){
    DEC_calc = PI*rad - Cart2Sph(U[0],U[1],U[2])[1]*rad;}


if (DEC_set<0){    
            DEC_calc =  Cart2Sph(U[0],U[1],U[2])[1]*rad - PI*rad;}


float Tot_er = abs(Angsepsph(RA_set,DEC_set,RA_calc/rad,DEC_calc/rad));

}

void Foy(float Pini[6],float Err, float RAs, float DECs, float RAest, float DECest){

    double Fin_Err_All[1000] = {};
    double Err_All[1000] = {};
    double RA_All[1000] = {};
    double DEC_All[1000] = {};

    if (DECs<0){
            DECs = 180 + DECs;}
    
    float Iter = {};

    float r = RAs;

    float D1[3] = {Pini[0]};
    float D2[3] = {Pini[1]};
    float D3[3] = {Pini[2]};
    float D4[3] = {Pini[3]};
    float D5[3] = {Pini[4]};
    float D6[3] = {Pini[5]};

    float P[6] = {D1[0],D2[0],D3[0],D4[0],D5[0],D6[0]};

    float P_first[3] = {D1[0],D1[1],D1[2]};
    
    D1[0] = D1[0] - P_first[0];
    D1[1] = D1[1] - P_first[1];
    D1[2] = D1[2] - P_first[2];

    D2[0] = D2[0] - P_first[0];
    D2[1] = D2[1] - P_first[1];
    D2[2] = D2[2] - P_first[2];

    D3[0] = D3[0] - P_first[0];
    D3[1] = D3[1] - P_first[1];
    D3[2] = D3[2] - P_first[2];

    D4[0] = D4[0] - P_first[0];
    D4[1] = D4[1] - P_first[1];
    D4[2] = D4[2] - P_first[2];

    D5[0] = D5[0] - P_first[0];
    D5[1] = D5[1] - P_first[1];
    D5[2] = D5[2] - P_first[2];

    D6[0] = D6[0] - P_first[0];
    D6[1] = D6[1] - P_first[1];
    D6[2] = D6[2] - P_first[2];

    float M = 6;

    RAs = RAs/rad;
    DECs = DECs/rad;

    float xs = r*cos(RAs)*sin(DECs);
    float ys = r*sin(RAs)*sin(DECs);
    float zs = r*cos(DECs);
    
    xs = xs -  P_first[0];
    ys = ys -  P_first[1];
    zs = zs -  P_first[2];

    float p_T[3] = {xs,ys,zs};

    float Ang[5] = {};
    
    Ang[0] = Angsepcart(D2[0],D2[1],D2[2],xs,ys,zs)  ;
    Ang[1] = Angsepcart(D3[0],D3[1],D3[2],xs,ys,zs)  ;
    Ang[2] = Angsepcart(D4[0],D4[1],D4[2],xs,ys,zs)  ;
    Ang[3] = Angsepcart(D5[0],D5[1],D5[2],xs,ys,zs)  ;
    Ang[4] = Angsepcart(D6[0],D6[1],D6[2],xs,ys,zs)  ;

    float toa1[3] = {xs-D1[0],ys-D1[1],zs-D1[2]} ;
    float toa2[3] = {xs-D2[0],ys-D2[1],zs-D2[2]} ;
    float toa3[3] = {xs-D3[0],ys-D3[1],zs-D3[2]} ;
    float toa4[3] = {xs-D4[0],ys-D4[1],zs-D4[2]} ;
    float toa5[3] = {xs-D5[0],ys-D5[1],zs-D5[2]} ;
    float toa6[3] = {xs-D6[0],ys-D6[1],zs-D6[2]} ;

    float tdoa[5] = {};
    tdoa[0] = toa2 - toa1;
    tdoa[1] = toa3 - toa1;
    tdoa[2] = toa4 - toa1;
    tdoa[3] = toa5 - toa1;
    tdoa[4] = toa6 - toa1;

    float xar = xs;
    float yar = ys;
    float zar = zs;

    float d_all[6]  = {};
    d_all[0] = length(xs,ys,zs,D1[0],D1[1],D1[2]);
    d_all[1] = length(xs,ys,zs,D2[0],D2[1],D2[2]);
    d_all[2] = length(xs,ys,zs,D3[0],D3[1],D3[2]);
    d_all[3] = length(xs,ys,zs,D4[0],D4[1],D4[2]);
    d_all[4] = length(xs,ys,zs,D5[0],D5[1],D5[2]);
    d_all[5] = length(xs,ys,zs,D6[0],D6[1],D6[2]);

    float t_all[5] = {};
    float c =299792458;
    t_all[0] = (d_all[1] - d_all[0])/c;
    t_all[1] = (d_all[2] - d_all[0])/c;
    t_all[2] = (d_all[3] - d_all[0])/c;
    t_all[3] = (d_all[4] - d_all[0])/c;
    t_all[4] = (d_all[5] - d_all[0])/c;
    
    /*  After Time Delay          */
    
    // Define Initial Guess

    float r_est = 10;
    float Vals[3] = {};
    float A = Sph2Cart(r_est,DECest/rad,RAest/rad,Vals);
    float p_T_0[3] = {Vals[0],Vals[1],Vals[2]};
    float Initial_guess[3] = {Vals[0],Vals[1],Vals[2]};

    // Add Uncertainty to Time-delay
    float TDOA[5] = {t_all[0],t_all[1],t_all[2],t_all[3],t_all[4]};

    float Err_0[5] = {Err,Err,Err,Err,0};

    int len[5] = {0, 1, 2, 3, 4};
    for( int y : len ) {
        TDOA[y] = TDOA[y] + Err_0[y];

    }

    double d[5] = {};
     for( int y : len ) {
        d[y] = TDOA[y]*c;

    }

   
    double i = 1;
    double f[5] = {0,0,0,0,0};
    double dDEC_f[5][3] = {};
    
    std::vector<int> v = { 1, 2, 3, 4, 5};

    double P_2[6][3] = {{D1[0],D1[1],D1[2]},{D2[0],D2[1],D2[2]},{D3[0],D3[1],D3[2]},{D4[0],D4[1],D4[2]},{D5[0],D5[1],D5[2]},{D6[0],D6[1],D6[2]}};
    
    // Convert to matrices
    double p_T_0_2[3] = {{p_T_0[0]},{p_T_0[1]},{p_T_0[2]}};

    for ( auto i : { 1,2,3,4,5,6 } )
        {
    
        double result1[3] = {};
        subtract2( p_T_0_2,P_2[i]  ,result1);
        float norm1 = norm(result1[0],result1[1],result1[2])    ;     
        
        double result2[3] = {};
        subtract2( p_T_0_2,P_2[0]  ,result2);
        float norm2 = norm(result2[0],result2[1],result2[2])    ;    

        f[i-1] = norm1 - norm2;

        dDEC_f[i-1][0] = (p_T_0_2[0] -P_2[i][0]) * pow(norm1,-1) - 
        (p_T_0_2[0] -P_2[0][0]) * pow(norm2,-1);


        dDEC_f[i-1][1] = (p_T_0_2[1] -P_2[i][1]) * pow(norm1,-1) - 
        (p_T_0_2[1] -P_2[0][1]) * pow(norm2,-1);

        dDEC_f[i-1][2] = (p_T_0_2[2] -P_2[i][2]) * pow(norm1,-1) - 
        (p_T_0_2[2] -P_2[0][2]) * pow(norm2,-1);

        
        }

    matrix dDEC_f2 = {{dDEC_f[0][0],dDEC_f[0][1],dDEC_f[0][2]},{dDEC_f[1][0],dDEC_f[1][1],dDEC_f[1][2]},{dDEC_f[2][0],dDEC_f[2][1],dDEC_f[2][2]},{dDEC_f[3][0],dDEC_f[3][1],dDEC_f[3][2]},{dDEC_f[4][0],dDEC_f[4][1],dDEC_f[4][2]}};

    matrix dDEC_f2inv = inverse(dDEC_f2);

    std::vector<double> input123({ dDEC_f2[0] });
    double dDEC_f2inv1[1];
    std::copy(input123.begin(), input123.end(), dDEC_f2inv1);

    std::vector<double> input124({ dDEC_f2[1] });
    double dDEC_f2inv2[1];
    std::copy(input124.begin(), input124.end(), dDEC_f2inv2);

    std::vector<double> input125({ dDEC_f2[2] });
    double dDEC_f2inv3[1];
    std::copy(input125.begin(), input125.end(), dDEC_f2inv3);

    double dDEC_f2invfin[3] = {dDEC_f2inv1[0],dDEC_f2inv2[1],dDEC_f2inv3[2]};

    double result3[3] = {};
    subtract2( d,f  ,result3);
       
    double x_nonlin[3] = {};
    double fin[3] = {dot_product2(result3, dDEC_f2invfin    ) } ;
    
    add2(fin,p_T_0_2,x_nonlin)  ;

    float RA = Cart2Sph(x_nonlin[0],x_nonlin[1],x_nonlin[2])[2];
    float DEC = Cart2Sph(x_nonlin[0],x_nonlin[1],x_nonlin[2])[1];
    
    if (DECs< 90/rad){
            DEC = DEC - 90/rad;}
    else{
        DEC = 270/rad - DEC ;}
       
    double Fin_Err = abs(Angsepsph(RAs,DECs,RA,DEC))    ;
    Fin_Err_All[0] = (Fin_Err);

    for ( auto i : { 1,2,3,4,5,6 } ){
        double Initial_guess[3] = {x_nonlin[0],x_nonlin[1],x_nonlin[2]};
        double p_T_0_3[3] = {x_nonlin[0],x_nonlin[1],x_nonlin[2]};

        double f[5] = {0,0,0,0,0};
        double dDEC_f[5][3] = {};
    
        double p_T_0_2[3] = {{x_nonlin[0]},{x_nonlin[1]},{x_nonlin[2]}};

        for ( auto i : { 1,2,3,4,5,6 } )
            {

            double result1[3] = {};
            subtract2( p_T_0_2,P_2[i]  ,result1);
            float norm1 = norm(result1[0],result1[1],result1[2])    ;     
            
            double result2[3] = {};
            subtract2( p_T_0_2,P_2[0]  ,result2);
            float norm2 = norm(result2[0],result2[1],result2[2])    ;    

            f[i-1] = norm1 - norm2;

            dDEC_f[i-1][0] = (p_T_0_2[0] -P_2[i][0]) * pow(norm1,-1) - 
            (p_T_0_2[0] -P_2[0][0]) * pow(norm2,-1);


            dDEC_f[i-1][1] = (p_T_0_2[1] -P_2[i][1]) * pow(norm1,-1) - 
            (p_T_0_2[1] -P_2[0][1]) * pow(norm2,-1);

            dDEC_f[i-1][2] = (p_T_0_2[2] -P_2[i][2]) * pow(norm1,-1) - 
            (p_T_0_2[2] -P_2[0][2]) * pow(norm2,-1);

            
            }

        matrix dDEC_f2 = {{dDEC_f[0][0],dDEC_f[0][1],dDEC_f[0][2]},{dDEC_f[1][0],dDEC_f[1][1],dDEC_f[1][2]},{dDEC_f[2][0],dDEC_f[2][1],dDEC_f[2][2]},{dDEC_f[3][0],dDEC_f[3][1],dDEC_f[3][2]},{dDEC_f[4][0],dDEC_f[4][1],dDEC_f[4][2]}};

        matrix dDEC_f2inv = inverse(dDEC_f2);

        std::vector<double> input123({ dDEC_f2[0] });
        double dDEC_f2inv1[1];
        std::copy(input123.begin(), input123.end(), dDEC_f2inv1);

        std::vector<double> input124({ dDEC_f2[1] });
        double dDEC_f2inv2[1];
        std::copy(input124.begin(), input124.end(), dDEC_f2inv2);

        std::vector<double> input125({ dDEC_f2[2] });
        double dDEC_f2inv3[1];
        std::copy(input125.begin(), input125.end(), dDEC_f2inv3);


        double dDEC_f2invfin[3] = {dDEC_f2inv1[0],dDEC_f2inv2[1],dDEC_f2inv3[2]};

        double result3[3] = {};
        subtract2( d,f  ,result3);
           
        double x_nonlin[3] = {};
        double fin[3] = {dot_product2(result3, dDEC_f2invfin    ) } ;
        
        add2(fin,p_T_0_2,x_nonlin)  ;

        float RA = Cart2Sph(x_nonlin[0],x_nonlin[1],x_nonlin[2])[2];
        float DEC = Cart2Sph(x_nonlin[0],x_nonlin[1],x_nonlin[2])[1];
        
        if (DECs< 90/rad){
                DEC = DEC - 90/rad;}
        else{
            DEC = 270/rad - DEC ;}
        
        
        double Fin_Err = abs(Angsepsph(RAs,DECs,RA,DEC))    ;
        Fin_Err_All[i] = (Fin_Err);

        Err_All[i] = Fin_Err;
        RA_All[i] = RA;
        DEC_All[i] = DEC;

    }
    


}

int main(){}





