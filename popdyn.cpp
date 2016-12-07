#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include "./Random_Number_Generators/randgen.h"

using namespace std;
using namespace randgenbase;

#include "./Random_Number_Generators/randgen.cpp"

int find_descriptor(ifstream& fp_in, const string& descriptor){
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
int probrange_check(const double& value, const string name){
	if((value <= double(0.0)) || (value > double(1.0)){
		cout << "Input error, value out of range for double variable " <<  \
		  name << endl;
		return(1);
	}
	return(0);
}
int positive_longcheck(const long& value, const string name){
	if(value <= 0){
		cout << "Input error, value out of range for integer variable " << \
		  name << endl;
		return(1);
	}
	return(0);
}
long direct_binomial(const double& success_prob, const long& trials){
	long count(0);
	for(long i=0; i<trials; ++i){
		if(randgen() < success_prob){
			++count;
		}
	}
	return(count);
}
long gauss_binomial(const double& success_prob, const long& trials){
	const double mean_success = static_cast<double>(trials)*success_prob;
	const double stddev = sqrt(mean_success*(double(1.0)-success_prob));
	const count = static_cast<long>(round(stddev*gaussrandgen() + \
	  mean_success));
	return(count);
}
long poisson_binomial(const double& success_prob, const long& trials){
	return(0);
}


int main(int argc, char** argv){
	if(argc < 2){
		cout << "No input file specified." << endl;
		exit(1);
	}
	string inputfile = argv[1];
	ifstream fp_in;
	int input_success(0);
	double mu(0.0),sigma(0.0),lambda(0.0);
	long max_timesteps(0),sample_interval(0);
	long prey_start_number(0),predator_start_number(0);
	
	fp_in.open(inputfile, ios::in);
	fp_in.clear();
	fp_in.seekg(0, ios::beg);
	input_success += find_descriptor(fp_in, "max_timesteps=");
	fp_in >> max_timesteps;
	input_success += positive_longcheck(max_timesteps, "max_timesteps");
	input_success += find_descriptor(fp_in, "sample_interval=");
	fp_in >> sample_interval;
	input_success += positive_longcheck(sample_interval, "sample_interval");
	input_success += find_descriptor(fp_in, "mu=");
	fp_in >> mu;
	input_success += probrange_check(mu, "mu");
	input_success += find_descriptor(fp_in, "sigma=");
	fp_in >> sigma;
	input_success += probrange_check(sigma, "sigma");
	input_success += find_descriptor(fp_in, "lambda=");
	fp_in >> lambda;
	input_success += probrange_check(lambda, "lambda");
	input_success += find_descriptor(fp_in, "prey_start_number=");
	fp_in >> prey_start_number;
	input_success += positive_longcheck(prey_start_number, \
	  "prey_start_number");
	input_success += find_descriptor(fp_in, "predator_start_number=");
	fp_in >> predator_start_number;
	input_success += positive_longcheck(predator_start_number, \
	  "predator_start_number");
	fp_in.close();
	if(input_success != 0){
		cout << "Error, required input parameters  ";
		cout << "missing or out of bounds." << endl;
		return(1);
	}
	
	const long max_samples = static_cast<long>( \
	  floor(static_cast<double>(max_timesteps)/ \
	  static_cast<double>(sample_interval)))+2;
	const int species = 2;
	long** samples = new long*[max_samples];
	for(int i=0; i<max_samples; ++i){
		samples[i] = new long[species];
		for(int j=0; j<species; ++j){
			samples[i][j] = 0;
		}
	}
	long* species_counts = new long[species];
	for(int i=0; i<species; ++i){
		species_counts[i] = 0;
	}
	autoinit_randgen();
	//const long max_relevant_timesteps = max_samples*sample_interval;
	const long max_relevant_timesteps = max_timesteps;
	const long fixed_sample_interval = sample_interval;
	long last_timestep(0),trigger_index(0),sample_index(1);
	species_counts[0] = prey_start_number;
	species_counts[1] = predator_start_number;
	samples[0][0] = species_counts[0];
	samples[0][1] = species_counts[1];
	for(long i=0; i<max_relevant_timesteps; ++i){
		long mu_change = direct_binomial(mu, species_counts[0]);
		long sigma_change = direct_binomial(sigma, species_counts[1]);
		long species_product = species_counts[0]*species_counts[1];
		long lambda_change = direct_binomial(lambda, species_product);
		species_counts[0] += (mu_change - lambda_change);
		species_counts[1] += (lambda_change - sigma_change);
		if((species_counts[0] <= 0) || (species_counts[1] <= 0)){
			last_timestep = i+1;
			break;
		}
		++trigger_index;
		if(trigger_index == fixed_sample_interval){
			samples[sample_index][0] = species_counts[0];
			samples[sample_index][1] = species_counts[1];
			++sample_index;
			trigger_index = 0;
		}
	}
	cout << endl;
	if(last_timestep > 0){
		cout << "Population extinction at timestep: " << last_timestep << endl;
	}
	else{
		cout << "Simulation finished with populations still existing." << endl;
		cout << "Simulated timesteps: " << max_relevant_timesteps << endl;
	}
	for(int i=0; i<species; ++i){
		if(species_counts[i] < 0){
			species_counts[i] = 0;
		}
	}
	cout << "Final populations: " << species_counts[0] << \
	  "  " << species_counts[1] << endl;
	cout << endl;
	ofstream fp_out;
	fp_out.open("population_trajectory.txt", ios::out);
	fp_out.clear();
	fp_out.seekp(0, ios::beg);
	for(long i=0; i<sample_index; ++i){
		fp_out << i*fixed_sample_interval << " " << samples[i][0] << " " << samples[i][1] << endl;
	}
	if(last_timestep > 0){
		fp_out << last_timestep << " " << species_counts[0] << " " << species_counts[1] << endl;
	}
	fp_out.close();
	
	if(species_counts){
		delete[] species_counts;
		species_counts = 0;
	}
	if(samples){
		for(int i=0; i<max_samples; ++i){
			delete[] samples[i];
		}
		delete[] samples;
		samples = 0;
	}
	return(0);
}
