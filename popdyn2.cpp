#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include "./Random_Number_Generators/randgen.h"

using namespace std;
using namespace randgenbase;

#include "./Random_Number_Generators/randgen.cpp"

constexpr int species = 2;
constexpr int two_int = 2;

//The following block of functions assist in reading the input
//file and checking that input values make sense.
int find_descriptor(ifstream& fp_in, const string& descriptor) {
	//This function searches an input file for a keyword descriptor
	//which precedes an input value.
	std::string in_string = "empty";
	while (!fp_in.eof()) {
		fp_in >> in_string;
		if (in_string == descriptor) {
			return(0);
		}
	}
	fp_in.clear();
	fp_in.seekg(0, ios::beg);
	while (!fp_in.eof()) {
		fp_in >> in_string;
		if (in_string == descriptor) {
			return(0);
		}
	}
	cout << "find_descriptor, unable to find descriptor " \
	  << descriptor << endl;
	return(1);
}
int probrange_check(const double& value, const string& name) {
	//This function checks that a value is between 0 and 1,
	//as a probability value should be.
	if ((value <= 0.0) || (value > 1.0)) {
		cout << "Input error, value out of range for double variable " <<  \
		  name << endl;
		return(1);
	}
	return(0);
}
int positive_longcheck(const long& value, const string& name) {
	//This function checks that a long integer is positive.
	if (value <= 0) {
		cout << "Input error, value out of range for integer variable " << \
		  name << endl;
		return(1);
	}
	return(0);
}
int binary_intcheck(const int value, const string& name) {
	//This function checks that an integer is either 0 or 1 or 2
	//and not any other value.
	if ((value != 0) && (value != 1) && (value !=2)) {
		cout << "Input error, non-binary choice for " << name << endl;
		return(1);
	}
	return(0);
}

//The following block of functions assist with the simulation operation.
inline long direct_binomial(const double& success_prob, const long& trials) {
	//This function calculates a long integer value from 
	//a binomial distribution using the direct method.
	long count(0);
	for (long i=0; i<trials; ++i) {
		if (randgen() < success_prob) {
			++count;
		}
	}
	return(count);
}
inline void advance_species(long*const species_counts_in, \
  const double& mu_in, const double& sigma_in, const double& lambda_in){
	//This function advances the predator and prey counts, using the
	//input probabilities and the direct binomial probability calculation.
	const long mu_change = direct_binomial(mu_in, species_counts_in[0]);
	const long sigma_change = direct_binomial(sigma_in, species_counts_in[1]);
	const long species_product = species_counts_in[0]*species_counts_in[1];
	const long lambda_change = direct_binomial(lambda_in, species_product);
	species_counts_in[0] += (mu_change - lambda_change);				
	species_counts_in[1] += (lambda_change - sigma_change);
	return;
}
inline int determine_quadrant(const long*const species_counts_in, \
  const long& equil_prey_in, const long& equil_predator_in) {
	//This function determines what "quadrant" the
	//predator prey system is in.
	int quadrant_out(-1);
	if (species_counts_in[1] > equil_predator_in) {
		if (species_counts_in[0] > equil_prey_in) {
			quadrant_out = 2;
		}
		else {
			quadrant_out = 1;
		}
	}
	else{
		if (species_counts_in[0] > equil_prey_in) {
			quadrant_out = 3;
		}
		else{
			quadrant_out = 0;
		}
	}
	return(quadrant_out);
}


int main(int argc, char** argv) {
	if(argc < 2){
		cout << "No input file specified." << endl;
		exit(1);
	}

	//The following block creates variables for reading in the
	//input file parameters.
	string inputfile = argv[1];
	ifstream fp_in;
	double mu_input(0.0);
	double sigma_input(0.0);
	double lambda_input(0.0);
	long max_timesteps_input(0);
	long sample_interval_input(0);
	long prey_start_number_input(0);
	long predator_start_number_input(0);
	int gather_extinction_times(0);
	string extinction_time_data = "empty";
	int keep_all_extinction_times(0);
	int keep_extinction_time_histogram(0);
	long simulation_trials(0);
	long extinction_time_histogram_length(0);
	long extinction_time_histogram_bin_size(1);
	int input_success(0);
	
	//The following block reads the input file and checks
	//the input parameters.
	fp_in.open("inputfile.txt", ios::in);
	fp_in.clear();
	fp_in.seekg(0, ios::beg);
	input_success += find_descriptor(fp_in, "max_timesteps=");
	fp_in >> max_timesteps_input;
	input_success += positive_longcheck(max_timesteps_input, "max_timesteps");
	input_success += find_descriptor(fp_in, "sample_interval=");
	fp_in >> sample_interval_input;
	input_success += positive_longcheck(sample_interval_input, \
	  "sample_interval");
	input_success += find_descriptor(fp_in, "mu=");
	fp_in >> mu_input;
	input_success += probrange_check(mu_input, "mu");
	input_success += find_descriptor(fp_in, "sigma=");
	fp_in >> sigma_input;
	input_success += probrange_check(sigma_input, "sigma");
	input_success += find_descriptor(fp_in, "lambda=");
	fp_in >> lambda_input;
	input_success += probrange_check(lambda_input, "lambda");
	input_success += find_descriptor(fp_in, "prey_start_number=");
	fp_in >> prey_start_number_input;
	if (prey_start_number_input != -1) {
		input_success += positive_longcheck(prey_start_number_input, \
		  "prey_start_number");
	}
	input_success += find_descriptor(fp_in, "predator_start_number=");
	fp_in >> predator_start_number_input;
	if (predator_start_number_input != -1) {
		input_success += positive_longcheck(predator_start_number_input, \
		  "predator_start_number");
	}
	input_success += find_descriptor(fp_in, "gather_extinction_times=");
	fp_in >> gather_extinction_times;
	input_success += binary_intcheck(gather_extinction_times, \
	  "gather_extinction_times");
	if (gather_extinction_times == 1) {
		input_success += find_descriptor(fp_in, "simulation_trials=");
		fp_in >> simulation_trials;
		input_success += positive_longcheck(simulation_trials, \
		  "simulation_trials");
		input_success += find_descriptor(fp_in, "extinction_time_data=");
		fp_in >> extinction_time_data;
		if (extinction_time_data == "times") {
			keep_all_extinction_times = 1;
		}
		else if (extinction_time_data == "histogram") {
			keep_extinction_time_histogram = 1;
		}
		else if (extinction_time_data == "all") {
			keep_all_extinction_times = 1;
			keep_extinction_time_histogram = 1;
		}
		else {
			cout << "Unknown extinction_time_data value: " << \
			  extinction_time_data << endl;
			++input_success;
		}
		if (keep_extinction_time_histogram == 1) {
			input_success += find_descriptor(fp_in, \
			  "extinction_time_histogram_length=");
			fp_in >> extinction_time_histogram_length;
			input_success += positive_longcheck(extinction_time_histogram_length, \
			  "extinction_time_histogram_length");
			input_success += find_descriptor(fp_in, \
			  "extinction_time_histogram_bin_size=");
			fp_in >> extinction_time_histogram_bin_size;
			input_success += positive_longcheck( \
			  extinction_time_histogram_bin_size, \
			  "extinction_time_histogram_bin_size");
		}
	}
	fp_in.close();
	if (input_success != 0) {
		cout << "Error, required input parameters  ";
		cout << "missing or out of bounds." << endl;
		return(1);
	}
	
	//The following block prepares the data structures that
	//will hold the gathered simulation data.
	const long max_samples = static_cast<long>( \
	  floor(static_cast<double>(max_timesteps_input)/ \
	  static_cast<double>(sample_interval_input))) + 2;
	const long max_sample_entries = max_samples*species;
	long* extinction_times = nullptr;
	long* extinction_time_histogram = nullptr;
	long* max_cycles_list = nullptr;
	long* samples = nullptr;
	if (gather_extinction_times == 0) {
		samples = new long[max_sample_entries];
		fill(samples, samples + max_sample_entries, -1);
	}
	else if (gather_extinction_times == 1) {
		max_cycles_list = new long[simulation_trials];
		fill(max_cycles_list, max_cycles_list + simulation_trials, -1);
		if (keep_all_extinction_times == 1) {
			extinction_times = new long[simulation_trials];
			fill(extinction_times, extinction_times + simulation_trials, -1);
		}
		if (keep_extinction_time_histogram == 1) {
			extinction_time_histogram = \
			  new long[extinction_time_histogram_length];
			fill(extinction_time_histogram, extinction_time_histogram + \
			  extinction_time_histogram_length, 0);
		}
	}

	//This block uses the input variables to prepare simulation variables.
	//It creates const variables for anything that will not change
	//during the simulation.
	const double mu = mu_input;
	const double sigma = sigma_input;
	const double lambda = lambda_input;
	const long equil_predator = static_cast<long>(floor(mu/lambda));
	const long equil_prey = static_cast<long>(floor(sigma/lambda));
	const long max_timesteps = max_timesteps_input;
	const long sample_interval = sample_interval_input;
	if (prey_start_number_input == -1) {
		prey_start_number_input = equil_prey;
		if (prey_start_number_input < 1) {
			prey_start_number_input = 1;
		}
	}
	if (predator_start_number_input == -1) {
		predator_start_number_input = equil_predator;
		if (predator_start_number_input < 1) {
			predator_start_number_input = 1;
		}
	}
	const long prey_start_number = prey_start_number_input;
	const long predator_start_number = predator_start_number_input;

	//This is the small array that holds the predator and prey counts.
	//It is updated as the numbers of predator and prey change.
	long species_counts[2] = {prey_start_number, predator_start_number};

	//This initializes the random number generator.
	autoinit_randgen();

	long cycle_count(0);
	long last_timestep(0);
	long quadrant = determine_quadrant(species_counts, equil_prey, \
	  equil_predator);
	long previous_quadrant = quadrant;
	const long start_quadrant = quadrant;
	if (gather_extinction_times == 0) {
		long trigger_index(0);
		long sample_index(0);
		samples[0] = species_counts[0];
		samples[1] = species_counts[1];
		sample_index = 1;
		for (long i=0; i<max_timesteps; ++i) {
			advance_species(species_counts, mu, sigma, lambda);
			for (int j=0; j<species; ++j) {
				if (species_counts[j] < 1) {
					species_counts[j] = 0;
					last_timestep = i + 1;
				}
			}
			quadrant = determine_quadrant(species_counts, equil_prey, \
			  equil_predator);
			if ((previous_quadrant == 3) && (quadrant == 0)) {
				--cycle_count;
			}
			else if ((previous_quadrant == 0) && (quadrant == 3)){
				++cycle_count;
			}
			/*
			else if (abs(previous_quadrant - quadrant) == two_int) {
				cout << "quadrant jump, previous = " << previous_quadrant << \
				  ",  current = " << quadrant <<  ",  timestep = " << \
				  i << endl;
			}
			*/
			if (last_timestep != 0) {
				break;
			}
			previous_quadrant = quadrant;
			++trigger_index;
			if (trigger_index == sample_interval) {
				samples[species*sample_index] = species_counts[0];
				samples[(species*sample_index) + 1] = species_counts[1];
				++sample_index;
				trigger_index = 0;
			}
		}
		cout << endl;
		if (last_timestep > 0) {
			cout << "Population extinction at timestep: " << last_timestep \
			  << endl;
		}
		else {
			cout << "Simulation finished with populations still existing." \
			  << endl;
			cout << "Simulated timesteps: " << max_timesteps << endl;
		}
		cout << "Final populations: " << endl;
		cout << "Prey: " << species_counts[0] << endl; 
		cout << "Predators: " << species_counts[1] << endl;
		cout << "Total population cycles: " << cycle_count << endl;
		cout << endl;
		ofstream fp_out;
		fp_out.open("population_trajectory.txt", ios::out);
		fp_out.clear();
		fp_out.seekp(0, ios::beg);
		fp_out << "timestep,prey,predators" << endl;
		for (long i=0; i<sample_index; ++i) {
			fp_out << i*sample_interval << "," << \
			  samples[i*species] << "," << \
			  samples[(i*species) + 1] << endl;
		}
		if (last_timestep > 0) {
			fp_out << last_timestep << "," << species_counts[0] << ","  \
			  << species_counts[1] << endl;
		}
		fp_out.close();
	}
	else if (gather_extinction_times == 1) {
		const double hbin_size = \
		  static_cast<double>(extinction_time_histogram_bin_size);
		for (long trial=0; trial<simulation_trials; ++trial) {
			cout << "trial " << trial << endl;
			cycle_count = 0;
			last_timestep = 0;
			species_counts[0] = prey_start_number;
			species_counts[1] = predator_start_number;
			quadrant = start_quadrant;
			previous_quadrant = start_quadrant;
			for (long i=0; i<max_timesteps; ++i) {
				advance_species(species_counts, mu, sigma, lambda);
				for (int j=0; j<species; ++j) {
					if (species_counts[j] < 1) {
						species_counts[j] = 0;
						last_timestep = i + 1;
					}
				}
				quadrant = determine_quadrant(species_counts, equil_prey, \
				  equil_predator);
				if ((previous_quadrant == 3) && (quadrant == 0)) {
					--cycle_count;
				}
				else if ((previous_quadrant == 0) && (quadrant == 3)) {
					++cycle_count;
				}
				if (last_timestep != 0) {
					break;
				}
				previous_quadrant = quadrant;
			}
			if (last_timestep == 0) {
				cout << "Warning, no extinction before timestep limit." << endl;
			}
			else {
				/*cout << "cycle count = " << cycle_count << endl;*/
				max_cycles_list[trial] = cycle_count;
				if (keep_all_extinction_times == 1) {
					extinction_times[trial] = last_timestep;
					/*cout << "last timestep = " << last_timestep << endl;*/
				}
				if (keep_extinction_time_histogram == 1) {
					long index = static_cast<long>( \
					round(static_cast<double>(last_timestep)/hbin_size));
					if(index < extinction_time_histogram_length){
						++extinction_time_histogram[index];
					}
					else {
						cout << "Warning, extinction time exceeds ";
						cout << "histogram length." << endl;
					}
				}
			}
		}
    }
	else if (gather_extinction_times == 2) {
			const double hbin_size = \
		  static_cast<double>(extinction_time_histogram_bin_size);
          long timestep_average = 0;
          long double timestep_sum = 0;
          long trial_count = 0;
		for (long trial=0; trial<simulation_trials; ++trial) {
            ++trial_count;
			cycle_count = 0;
			last_timestep = 0;
			species_counts[0] = prey_start_number;
			species_counts[1] = predator_start_number;
			quadrant = start_quadrant;
			previous_quadrant = start_quadrant;
			for (long i=0; i<max_timesteps; ++i) {
				advance_species(species_counts, mu, sigma, lambda);
				for (int j=0; j<species; ++j) {
					if (species_counts[j] < 1) {
						species_counts[j] = 0;
						last_timestep = i + 1;
					}
				}
				quadrant = determine_quadrant(species_counts, equil_prey, \
				  equil_predator);
				if ((previous_quadrant == 3) && (quadrant == 0)) {
					--cycle_count;
				}
				else if ((previous_quadrant == 0) && (quadrant == 3)) {
					++cycle_count;
				}
				if (last_timestep != 0) {
					break;
				}
				previous_quadrant = quadrant;
			}
			if (last_timestep == 0) {
				cout << "Warning, no extinction before timestep limit." << endl;
			}
			else {
				/*cout << "cycle count = " << cycle_count << endl;*/
				max_cycles_list[trial] = cycle_count;
				if (keep_all_extinction_times == 1) {
					extinction_times[trial] = last_timestep;
					/*cout << "last timestep = " << last_timestep << endl;*/
				}
				if (keep_extinction_time_histogram == 1) {
					long index = static_cast<long>( \
					round(static_cast<double>(last_timestep)/hbin_size));
					if(index < extinction_time_histogram_length){
						++extinction_time_histogram[index];
					}
					else {
						cout << "Warning, extinction time exceeds ";
						cout << "histogram length." << endl;
					}
				}
            }
        }
        timestep_sum = timestep_sum + last_timestep;
        timestep_average = timestep_sum / trial_count;
        cout << "Final populations: " << endl;
		cout << "Prey: " << species_counts[0] << endl; 
		cout << "Predators: " << species_counts[1] << endl;
		cout << "Total population cycles: " << cycle_count << endl;
		cout << endl;
        cout << "Average timestep to extinction: " << timestep_average << endl;
        cout << endl;
    }
		ofstream fp_out;
		fp_out.open("max_cycles_list.txt", ios::out);
		fp_out.clear();
		fp_out.seekp(0, ios::beg);
		for (long i=0; i<simulation_trials; ++i) {
			if (max_cycles_list[i] > -1) {
				fp_out << max_cycles_list[i] << endl;
			}
		}
		fp_out.close();
		if (keep_all_extinction_times == 1) {
			fp_out.open("extinction_time_list.txt", ios::out);
			fp_out.clear();
			fp_out.seekp(0, ios::beg);
			for(long i=0; i<simulation_trials; ++i){
				if(extinction_times[i] > 0){
					fp_out << extinction_times[i] << endl;
				}
			}
			fp_out.close();
		}
		if (keep_extinction_time_histogram == 1) {
			fp_out.open("extinction_time_histogram.txt", ios::out);
			fp_out.clear();
			fp_out.seekp(0, ios::beg);
			for (long i=0; i<extinction_time_histogram_length; ++i) {
				fp_out << i*extinction_time_histogram_bin_size << " " \
				  << extinction_time_histogram[i] << endl;
			}
			fp_out.close();
		}
	    
	if (samples) {
		delete[] samples;
		samples = nullptr;
	}
	if (extinction_times) {
		delete[] extinction_times;
		extinction_times = nullptr;
	}
	if (extinction_time_histogram) {
		delete[] extinction_time_histogram;
		extinction_time_histogram = nullptr;
	}
	if (max_cycles_list) {
		delete[] max_cycles_list;
		max_cycles_list = nullptr;
	}
	return(0);
}
