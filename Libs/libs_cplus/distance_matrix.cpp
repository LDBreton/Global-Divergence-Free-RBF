#include <iostream>
#include "Eigen/MPRealSupport"
#include "Eigen/LU"
#include <omp.h>
#include <vector>

using namespace mpfr;
using namespace Eigen;
using namespace std;  

typedef Matrix<mpreal,Dynamic,Dynamic>  MatrixXmp;
//typedef mpreal (*FnPtr)(mpreal,mpreal,mpreal,mpreal);


template<typename TFunc>
MatrixXmp distance_matrix (const vector< vector<mpreal> >  Puntos, const vector< vector<mpreal> >  Centros, const TFunc fbr,const vector< mpreal > paramrs)
{
  int Npuntos = Puntos.size();
  int NCentros = Centros.size();
  MatrixXmp fbrmat(Npuntos,NCentros);
  
  
  if(Npuntos == 0 || NCentros ==0 ){
	  return fbrmat;
	  }
  
  for(int i=0; i < Npuntos; i++){
    for(int j=0; j < NCentros; j++){
	fbrmat(i,j)=fbr(Puntos[i][0],Puntos[i][1],Centros[j][0],Centros[j][1],paramrs);		
	}
  }
  
  return fbrmat;
}

template<typename TFunc>
MatrixXmp distance_matrix (const MatrixXmp Puntos,const MatrixXmp Centros,const TFunc fbr,const vector< mpreal > paramrs)
{
  int Npuntos = Puntos.rows();
  int NCentros = Centros.rows();
  MatrixXmp fbrmat(Npuntos,NCentros);
  
  
  if(Npuntos == 0 || NCentros ==0 ){
	  return fbrmat;
	  }
  
  for(int i=0; i < Npuntos; i++){
    for(int j=0; j < NCentros; j++){
	fbrmat(i,j)=fbr(Puntos(i,0),Puntos(i,1),Centros(j,0),Centros(j,1),paramrs);		
	}
  }
  
  return fbrmat;
}



template<typename TFunc>
MatrixXmp distance_matrix_symetric (const vector< vector<mpreal> >  Puntos,const vector< vector<mpreal> >  Centros,const TFunc fbr,const vector< mpreal > paramrs)
{
  int Npuntos = Puntos.size();
  int NCentros = Centros.size();
  MatrixXmp fbrmat(Npuntos,NCentros);
  
  
  if(Npuntos == 0 || NCentros ==0 ){
	  return fbrmat;
	  }
  
  
  for(int i=0; i < Npuntos; i++){
    for(int j=0; j <= i; j++){
	fbrmat(i,j)=fbr(Puntos[i][0],Puntos[i][1],Centros[j][0],Centros[j][1],paramrs);
	fbrmat(j,i)=fbrmat(i,j); 		
	}
  }
  
  return fbrmat;
}


template<typename TFunc>
MatrixXmp distance_matrix_symetric (const MatrixXmp Puntos,const MatrixXmp Centros, const TFunc fbr,const vector< mpreal > paramrs)
{
  int Npuntos = Puntos.rows();
  int NCentros = Centros.rows();
  MatrixXmp fbrmat(Npuntos,NCentros);
  
  
  if(Npuntos == 0 || NCentros ==0 ){
	  return fbrmat;
	  }
  
  for(int i=0; i < Npuntos; i++){
    for(int j=0; j <= i ; j++){
	fbrmat(i,j)=fbr(Puntos(i,0),Puntos(i,1),Centros(j,0),Centros(j,1),paramrs);	
	fbrmat(j,i)=fbrmat(i,j); 		
	
	}
  }
  
  return fbrmat;
}
