#include <fstream>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <string>
#include "Eigen/MPRealSupport"

using namespace mpfr;
using namespace Eigen;
using namespace std;  

typedef Matrix<mpreal,Dynamic,Dynamic>  MatrixXmp;

//para leer los indices generados por matlab
vector< vector<int> > read_indexfile(string filename){
	    ifstream ifile(filename, std::ios::in);
	    int num;
	    string   line;
        vector< vector<int> > index;
		   while (getline(ifile, line)) {
		   vector<int> subindex;
           std::stringstream linestream(line);
           while (linestream >> num) {
				subindex.push_back(num);
			}
			index.push_back(subindex);
    }

	return index;
	
	}

//para leer archivo de puntos a mppreal
vector< vector<mpreal> > read_mpreal(string filename){
	    //mpreal::set_default_prec(mpfr::digits2bits(digits));
	    ifstream ifile(filename, std::ios::in);
	    string num;
	    string   line;
        vector< vector<mpreal> > mpreal_array;
		   while (getline(ifile, line)) {
		   vector<mpreal> subindex;
           std::stringstream linestream(line);
           while (linestream >> num) {
				subindex.push_back(num);
			}
			mpreal_array.push_back(subindex);
    }
	return mpreal_array;
	}
//para leer el archivo de puntos a matrices
MatrixXmp file2puntos(string filename){
	    vector< vector<mpreal> > aux = read_mpreal(filename);
	    int Npuntos = aux.size();
		MatrixXmp Puntos(Npuntos,2);
		
		for(int i=0; i < Npuntos; i++ ){
			for(int j=0; j < 2; j++ ){
				Puntos(i,j) = aux[i][j];
				}
			}
		return Puntos;

	}	


//para leer archivo de puntos a mppreal
vector<mpreal>  read_distance(string filename){
	    //mpreal::set_default_prec(mpfr::digits2bits(digits));
	    ifstream ifile(filename, std::ios::in);
	    string num;
        vector<mpreal> mpreal_distance;
           while (ifile >> num) {
				mpreal_distance.push_back(num);
			}
    
	return mpreal_distance;
	}	
