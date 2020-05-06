#include <iostream>
#include "Eigen/MPRealSupport"
#include "Eigen/LU"
#include <functional>   // std::bind
#include <map>
#include <omp.h>
#include <boost/typeof/typeof.hpp>
#include "distance_matrix.cpp"

using namespace mpfr;
using namespace Eigen;
using namespace std;

  
typedef Matrix<mpreal,Dynamic,Dynamic>  MatrixXmp;



template<typename Matrixx>
void insert_block(Matrixx &Mbig,Matrixx &Mblock,const int p, const int q){
	 int Nblock_cols = Mblock.cols();
     int Nblock_rows = Mblock.rows();
     
	if(Nblock_cols != 0 || Nblock_rows !=0 ){
     Mbig.block(p,q,Nblock_rows,Nblock_cols) = Mblock;
	  }
	}
	
	

template<typename ArrrayFunctions>
MatrixXmp BuildMGram(const vector< MatrixXmp > Puntos,const vector< MatrixXmp > Centros, const ArrrayFunctions Afbr_func,const vector< mpreal > paramrs){

int Nrows =0;
int Ncols =0;

	for(unsigned int i=0; i < Puntos.size() ; i++) {
	Nrows += Puntos[i].rows();
	}	
	for(unsigned int i=0; i < Centros.size() ; i++) {
	Ncols += Centros[i].rows();
	}
	MatrixXmp fbrgram_M(Nrows,Ncols);
	
	int p = 0;
	for(unsigned int i=0; i < Puntos.size() ; i++) {
		int q = 0;
        for(unsigned int j=0; j < Centros.size(); j++) {
            MatrixXmp fbrgramXX_M = distance_matrix(Puntos[i],Centros[j],Afbr_func[j][i],paramrs);
			insert_block(fbrgram_M,fbrgramXX_M,p,q);
		    q += Centros[j].rows();		     
        }
    p += Puntos[i].rows();    
	}
	
	return fbrgram_M;
}

template<typename ArrrayFunctions>
MatrixXmp BuildMGram_s(const vector< MatrixXmp > Puntos,const vector< MatrixXmp > Centros, const ArrrayFunctions Afbr_func,const vector< mpreal > paramrs){

int Nrows =0;
int Ncols =0;

	for(unsigned int i=0; i < Puntos.size() ; i++) {
	Nrows += Puntos[i].rows();
	}	
	for(unsigned int i=0; i < Centros.size() ; i++) {
	Ncols += Centros[i].rows();
	}
	MatrixXmp fbrgram_M(Nrows,Ncols);
	
	int p = 0;
	for(unsigned int i=0; i < Puntos.size() ; i++) {
		int q = 0;
        for(unsigned int j=0; j < Centros.size(); j++) {
			
			if(p<=q){
			if(p==q){
            MatrixXmp fbrgramXX_M = distance_matrix_symetric(Puntos[i],Centros[j],Afbr_func[j][i],paramrs);
    			insert_block(fbrgram_M,fbrgramXX_M,p,q);

			}
			else{
            MatrixXmp fbrgramXX_M = distance_matrix(Puntos[i],Centros[j],Afbr_func[j][i],paramrs);	
           			insert_block(fbrgram_M,fbrgramXX_M,p,q);

			}
			}
			
		    q += Centros[j].rows();		     
        }
    p += Puntos[i].rows();    
	}
	
	for(int i=0; i < Nrows; i++){
    for(int j=0; j <= i ; j++){
	fbrgram_M(j,i)=fbrgram_M(i,j); 		
	
	}
  }
	
	
	return fbrgram_M;
}


template<typename ArrrayFunctions>
MatrixXmp BuildMGram_t(const MatrixXmp Punto,const vector< MatrixXmp > Centros, const ArrrayFunctions Afbr_func,const vector< mpreal > paramrs){

int Nrows =0;

	for(unsigned int i=0; i < Centros.size() ; i++) {
	Nrows += Centros[i].rows();
	}
	
	MatrixXmp fbrgram_M(Nrows,1);
	int q = 0;
    
    for(unsigned int j=0; j < Centros.size(); j++) {
            MatrixXmp fbrgramXX_M = distance_matrix(Punto,Centros[j],Afbr_func[j],paramrs).transpose();
			insert_block(fbrgram_M,fbrgramXX_M,q,0);
		    q += Centros[j].rows();		     
        }
	
	return fbrgram_M;
}
