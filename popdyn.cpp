#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include "./Random_Number_Generators/randgen.h"

using namespace std;
using namespace randgenbase;

#include "./Random_Number_Generators/randgen.cpp"

int find_descriptor(ifstream& fp_in, const string descriptor){
	std::string in_string = "empty";
	while(!fp_in.eof()){
		fp_in >> in_string;
		if(in_string == descriptor){
			return(0);
		}
	}
	fp_in.clear();
	fp_in.seekg(0, ios::beg);
	while(!fp_in.eof()){
		fp_in >> in_string;
		if(in_string == descriptor){
			return(0);
		}
	}
	cout << "find_descriptor, unable to find descriptor " \
	  << descriptor << endl;
	return(1);
}



int main(int argc, char** argv){
	
	
	return(0);
}
