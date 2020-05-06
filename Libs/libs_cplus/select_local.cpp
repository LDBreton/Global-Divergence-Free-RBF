#include <fstream>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <string>
#include <Eigen/MPRealSupport>

using namespace mpfr;
using namespace Eigen;
using namespace std;  

typedef Matrix<mpreal,Dynamic,Dynamic>  MatrixXmp;


//por optimizar quizas clase?
//template<typename TFunc>
MatrixXmp select_local (const MatrixXmp Puntos,const vector<int> indexs)
{
  
  if(indexs[0] == -1){
	  MatrixXmp nulll(0,0);
	  return nulll;
	  }
  
  int Nsubpuntos = indexs.size(); 
  MatrixXmp subPuntos(Nsubpuntos,2);
    
  for(int i=0; i < Nsubpuntos; i++){
	  subPuntos(i,0) = Puntos(indexs[i],0);
	  subPuntos(i,1) = Puntos(indexs[i],1);
   }
  
  
  return subPuntos;
}
