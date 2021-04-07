#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
using namespace std;

// define a new short name for the STL vectors and vectors of vectors
using vect_t = std::vector<double>;
using matr_t = std::vector<std::vector<double>>;


// matrix multiplication using *
matr_t operator * (const matr_t &A, const matr_t &B)
{
  if (A.at(0).size() != B.size())
    throw std::runtime_error("First matrix must have the same number of columns as the rows of the second one.");

  size_t INT = B.size();
  size_t ROW = A.size();
  size_t COL = B.at(0).size();
  matr_t C(ROW, vect_t(COL));

	for (std::size_t row=0; row!=ROW; ++row)
 		for (std::size_t col=0; col!=COL; ++col)
      for (std::size_t i=0; i!=INT; ++i)   // alternatively, use std::inner_product
        C[row][col] += A[row][i] * B[i][col];

	return C;
}


const long double pi = 3.141592653589793238462643383279502884197169399375105820974944592307816406286209;
//const long double e = 2.7182818284590452353602874713526624977572470936999595749669676277240766303535476;
const long double q = 1.602176620898e-19;
const double m = 1.6726219e-27; // proton mass
const double c = 299792458;
const double L = 100; // lungh segmento -->  10T
const double deltaE = 1e6;  // E fornita nelle cavità acceleranti (eV?)
const double s = 0.5;  // length of the trajectory inside the quadropole


vect_t vect_mult_mat(matr_t const &B, vect_t const &V);
vect_t times2(vect_t const &V);
vect_t tubeStep(vect_t const &V);
vect_t accelStep(vect_t const &V);
vect_t focusStep(vect_t const &V);
vect_t defocusStep(vect_t const &V);
vect_t totalStep(vect_t V);
void print_mat(matr_t const &A);


int main()  {
  cout<<"\nHello\n\n";
  FILE *f = fopen("synchrotron.csv", "w");
  if (f == NULL) return -1;
  fprintf(f,"x,y,px,py,P\n");

  // initial phase space vector
  vect_t V0 {-0.1, 0, -3e-22, 5e-22, 1e-18};

  // print start vector
  cout<<"Starting particle state:" <<endl;
  for (auto &v : V0)
    cout << v << endl;
  cout<<endl;

  // Run of N/4 giri, write the trajectory in the .csv file
  int N = 2;
  for (int i=0; i!=N; ++i)
  {
    fprintf(f, "%0.18lf,%0.18lf,%0.18lf,%0.18lf,%0.18lf\n", V0[0], V0[1], V0[2], V0[3], V0[4]);
    V0 = totalStep(V0);
  }

  // print resulting vector
  cout<<endl<<"Final particle phase space vector after "<<N<<" blocks:"<<endl;
  for (auto &v : V0)
    cout << v << endl;

  cout<<endl;
  return 0;
}


vect_t totalStep(vect_t V)
{
  unsigned int dim = V.size();
  vect_t V1(dim);
  V1 = tubeStep(defocusStep(focusStep(tubeStep(accelStep(V)))));
  //matr_t Mtot(dim, vect_t(dim));
  return(V1);
}

vect_t tubeStep(vect_t const &V)
{
  //double P = sqrt( pow(V[2],2) + pow(V[3],2) + pow(V[4],2) );
  matr_t Mo { {1, 0, L/V[4], 0, 0},
              {0, 1, 0, L/V[4], 0},
              {0, 0, 1, 0, 0},
              {0, 0, 0, 1, 0},
              {0, 0, 0, 0, 1} };

  return vect_mult_mat(Mo, V);
}

vect_t accelStep(vect_t const &V)
{
  double P = sqrt( pow(V[2],2) + pow(V[3],2) + pow(V[4],2) );
  std::cout<<"Ptot = "<<P<<std::endl;
  matr_t Ma { {1, 0, 0, 0, 0},
              {0, 1, 0, 0, 0},
              {0, 0, 1, 0, 0},
              {0, 0, 0, 1, 0},
              {0, 0, 0, 0, 1 + deltaE/(c*P)} };

  return vect_mult_mat(Ma, V);
}

vect_t focusStep(vect_t const &V)
{
  //double P = sqrt( pow(V[2],2) + pow(V[3],2) + pow(V[4],2) );
  double B = V[4] / (q*2*s/pi); // p_z / (q * Raggio curvatura)
  //double B = 5;
  std::cout<<B<<"\n";
  double g = B / 0.1; // B / largh -->  ~10T/0.1m
  double k = q*g / V[4];
  double sqrk = sqrt(k);
  std::cout<<sqrk*s<<"\n";
  matr_t Mf { {cos(sqrk*s), 0, sin(sqrk*s)/(m*c*sqrk), 0, 0},
              {0, cosh(sqrk*s), 0, sinh(sqrk*s)/(m*c*sqrk), 0},
              {-m*c*sqrk*sin(sqrk*s), 0, cos(sqrk*s), 0, 0},
              {0, m*c*sqrk*sinh(sqrk*s), 0, cosh(sqrk*s), 0},
              {0, 0, 0, 0, 1 } };

  return vect_mult_mat(Mf, V);
}

vect_t defocusStep(vect_t const &V)
{
  //double P = sqrt( pow(V[2],2) + pow(V[3],2) + pow(V[4],2) );
  double B = V[4] / (q*2*s/pi);  // p_z / (q * Raggio curvatura)
  //double B = 5;
  double g = B / 0.1; // B / largh -->  ~10T/0.1m
  double k = q*g / V[4];
  double sqrk = sqrt(k);
  matr_t Md { {0, cosh(sqrk*s), 0, sinh(sqrk*s)/(m*c*sqrk), 0},
              {cos(sqrk*s), 0, sin(sqrk*s)/(m*c*sqrk), 0, 0},
              {0, m*c*sqrk*sinh(sqrk*s), 0, cosh(sqrk*s), 0},
              {-m*c*sqrk*sin(sqrk*s), 0, cos(sqrk*s), 0, 0},
              {0, 0, 0, 0, 1 } };

  return vect_mult_mat(Md, V);
}


vect_t vect_mult_mat(matr_t const &B, vect_t const &V)
{
  // aggiungere controllo dimensionale
  size_t dim = V.size();
  vect_t V2(dim);
  std::vector<double>::const_iterator col;  // declares a const iterator i for a vector<double>
  std::vector<double>::const_iterator row;
	for (col=V2.begin(); col!=V2.end(); ++col)
	{
    int coln = col - V2.begin();
 		for (row=V2.begin(); row!=V2.end(); ++row) {
      int rown = row - V2.begin();
      V2[coln] += B[coln][rown] * V.at(rown);  // .at() is slower but will give errors instead of segfault
    }
	}
	return V2;
}


void print_mat(const matr_t &A)
{
  for (auto &row : A) {
    for (auto &a : row) std::cout<< a <<"  ";
    std::cout<<std::endl;
  }
}
