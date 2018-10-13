#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <assert.h>
#include "./Random_Number_Generators/randgen.h"

using namespace std;
using namespace randgenbase;

#include "./Random_Number_Generators/randgen.cpp"

constexpr int species_types = 2;
constexpr double zero_double = 0.0;
constexpr int zero_int = 0;
constexpr long zero_long = 0;
constexpr long one_long = 1;
constexpr long neg_one_long = -1;

typedef struct simulation_parameters {
	//This struct holds the simulation input parameters.
	public:
	double mu;
	double sigma;
	double lambda;
	long equil_prey;
	long equil_predator;
	long prey_start_number;
	long predator_start_number;
	long max_timesteps;
	long simulation_trials;
	int keep_trajectory_samples;
	long sample_interval;
	int perform_start_number_sweep;
	int species_to_sweep;
	long sweep_count;
	long sweep_change;
	int keep_cycles_to_extinction;
	int keep_extinction_time_histogram;
	long extinction_time_histogram_length;
	long extinction_time_histogram_bin_size;
} simulation_parameters;

typedef struct simulation_data {
	//This is a struct to hold the gathered simulation data.
	public:
	long* extinction_times;
	long* extinction_time_histogram;
	long* sweep_average_extinction_times;
	long* max_cycles_list;
	long* samples;
	long single_group_samples;
	~simulation_data() {
		if (extinction_times) {
			delete[] extinction_times;
			extinction_times = nullptr;
		}
		if (extinction_time_histogram) {
			delete[] extinction_time_histogram;
			extinction_time_histogram = nullptr;
		}
		if (sweep_average_extinction_times) {
			delete[] sweep_average_extinction_times;
			sweep_average_extinction_times = nullptr;
		}
		if (max_cycles_list) {
			delete[] max_cycles_list;
			max_cycles_list = nullptr;
		}
		if (samples) {
			delete[] samples;
			samples = nullptr;
		}
		return;
	}
	simulation_data() {
		extinction_times = nullptr;
		extinction_time_histogram = nullptr;
		sweep_average_extinction_times = nullptr;
		max_cycles_list = nullptr;
		samples = nullptr;
		single_group_samples = zero_long;
		return;
	}
} simulation_data;


//The following block of functions assist in reading the input
//file and checking that input values make sense.
int find_descriptor(ifstream& fp_in, const string& descriptor) noexcept {
	//This function searches an input file for a keyword descriptor
	//which precedes an input value.
	if (!fp_in.is_open()) {
		cout << "find_descrioptor, file not open for descriptor " <<
		  descriptor << endl;
		return(1);
	}
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

int probrange_check(const double& value, const string& name) noexcept {
	//This function checks that a value is between 0 and 1,
	//as a probability value should be.
	if ((value <= 0.0) || (value > 1.0)) {
		cout << "Input error, value out of range for double variable " <<  \
		  name << endl;
		return(1);
	}
	return(0);
}

int positive_longcheck(const long& value, const string& name) noexcept {
	//This function checks that a long integer is positive.
	if (value <= 0) {
		cout << "Input error, value out of range for integer variable " << \
		  name << endl;
		return(1);
	}
	return(0);
}

int binary_intcheck(const int value, const string& name) noexcept {
	//This function checks that an integer is either 0 or 1
	//and not any other value.
	if ((value != 0) && (value != 1)) {
		cout << "Input error, non-binary choice for " << name << endl;
		return(1);
	}
	return(0);
}

int species_string_check(const string& value, const string& name) noexcept {
	if ((value != "predator") && (value != "prey")) {
		cout << "Input error, invalid species string choice for " << \
		  name << endl;
		return(1);
	}
	return(0);
}


int open_input_file(const string& filename, ifstream& fp_in, \
  const string& description) noexcept {
	if (fp_in.is_open()) {
		cout << "Error, input file already open for " << description << endl;
		return(1);
	}
	fp_in.open(filename, ios::in);
	if (!fp_in.is_open()) {
		cout << "Error, input file failed to open for " << description << endl;
		return(1);
	}
	fp_in.clear();
	fp_in.seekg(0, ios::beg);
	return(0);
}

int close_input_file(ifstream& fp_in, const string& description) noexcept {
	if (!fp_in.is_open()) {
		cout << "Error, input file already closed for " << description << endl;
		return(1);
	}
	fp_in.close();
	if (fp_in.is_open()) {
		cout << "Error, input file failed to close for " << description << \
		  endl;
		return(1);
	}
	return(0);
}

int open_output_file(const string& filename, ofstream& fp_out, \
  const string& description) noexcept {
	if (fp_out.is_open()) {
		cout << "Error, output file already open for " << description << endl;
		return(1);
	}
	fp_out.open(filename, ios::out);
	if (!fp_out.is_open()) {
		cout << "Error, output file failed to open for " << description << endl;
		return(1);
	}
	fp_out.clear();
	fp_out.seekp(0, ios::beg);
	return(0);
}

int close_output_file(ofstream& fp_out, const string& description) noexcept {
	if (!fp_out.is_open()) {
		cout << "Error, output file already closed for " << description << endl;
		return(1);
	}
	fp_out.close();
	if (fp_out.is_open()) {
		cout << "Error, output file failed to close for " << description << \
		  endl;
		return(1);
	}
	return(0);
}


int read_input_parameters(const string& filename, ifstream& fp_in, \
  simulation_parameters& params) noexcept {
	int input_success = 0;
	input_success += open_input_file(filename, fp_in, "base input file");
	input_success += find_descriptor(fp_in, "mu=");
	fp_in >> params.mu;
	input_success += probrange_check(params.mu, "mu");
	input_success += find_descriptor(fp_in, "sigma=");
	fp_in >> params.sigma;
	input_success += probrange_check(params.sigma, "sigma");
	input_success += find_descriptor(fp_in, "lambda=");
	fp_in >> params.lambda;
	input_success += probrange_check(params.lambda, "lambda");
	params.equil_prey = static_cast<long>(floor(params.sigma/params.lambda));
	params.equil_predator = static_cast<long>(floor(params.mu/params.lambda));
	input_success += find_descriptor(fp_in, "prey_start_number=");
	fp_in >> params.prey_start_number;
	if (params.prey_start_number == -1) {
		params.prey_start_number = params.equil_prey;
		if (params.prey_start_number < 1) {
			params.prey_start_number = 1;
		}
	}
	else {
		input_success += positive_longcheck(params.prey_start_number, \
		  "prey_start_number");
	}
	input_success += find_descriptor(fp_in, "predator_start_number=");
	fp_in >> params.predator_start_number;
	if (params.predator_start_number == -1) {
		params.predator_start_number = params.equil_predator;
		if (params.predator_start_number < 1) {
			params.predator_start_number = 1;
		}
	}
	else {
		input_success += positive_longcheck(params.predator_start_number, \
		  "predator_start_number");
	}
	input_success += find_descriptor(fp_in, "max_timesteps=");
	fp_in >> params.max_timesteps;
	input_success += positive_longcheck(params.max_timesteps, "max_timesteps");
	input_success += find_descriptor(fp_in, "simulation_trials=");
	fp_in >> params.simulation_trials;
	input_success += positive_longcheck(params.simulation_trials, \
	  "simulation_trials");
	input_success += find_descriptor(fp_in, "perform_start_number_sweep=");
	fp_in >> params.perform_start_number_sweep;
	input_success += binary_intcheck(params.perform_start_number_sweep, \
	  "perform_start_number_sweep");
	if (params.perform_start_number_sweep == 1) {
		input_success += find_descriptor(fp_in, "species_to_sweep=");
		string species_to_sweep_input = "";
		fp_in >> species_to_sweep_input;
		input_success += species_string_check(species_to_sweep_input, \
		  "species_to_sweep");
		const int chosen_species = (species_to_sweep_input == "prey") ? 0 : 1;
		params.species_to_sweep = chosen_species;
		input_success += find_descriptor(fp_in, "sweep_count=");
		fp_in >> params.sweep_count;
		input_success += positive_longcheck(params.sweep_count, "sweep_count");
		input_success += find_descriptor(fp_in, "sweep_change=");
		fp_in >> params.sweep_change;
	}
	input_success += find_descriptor(fp_in, "keep_trajectory_samples=");
	fp_in >> params.keep_trajectory_samples;
	input_success += binary_intcheck(params.keep_trajectory_samples, \
	  "keep_trajectory_samples");
	if (params.keep_trajectory_samples == 1) {
		input_success += find_descriptor(fp_in, "sample_interval=");
		fp_in >> params.sample_interval;
		input_success += positive_longcheck(params.sample_interval, \
		  "sample_interval");
	}
	input_success += find_descriptor(fp_in, "keep_cycles_to_extinction=");
	fp_in >> params.keep_cycles_to_extinction;
	input_success += binary_intcheck(params.keep_cycles_to_extinction, \
	  "keep_cycles_to_extinction");
	input_success += find_descriptor(fp_in, "keep_extinction_time_histogram=");
	fp_in >> params.keep_extinction_time_histogram;
	input_success += binary_intcheck(params.keep_extinction_time_histogram, \
	  "keep_extinction_time_histogram");
	if (params.keep_extinction_time_histogram == 1) {
		input_success += find_descriptor(fp_in, \
		  "extinction_time_histogram_length=");
		fp_in >> params.extinction_time_histogram_length;
		input_success += \
		  positive_longcheck(params.extinction_time_histogram_length, \
		  "extinction_time_histogram_length");
		input_success += find_descriptor(fp_in, \
		  "extinction_time_histogram_bin_size=");
		fp_in >> params.extinction_time_histogram_bin_size;
		input_success += positive_longcheck( \
		  params.extinction_time_histogram_bin_size, \
		  "extinction_time_histogram_bin_size");
	}
	input_success += close_input_file(fp_in, "base input file");
	if (input_success != 0) {
		cout << "Error, required input parameters  ";
		cout << "missing or out of bounds." << endl;
		return(1);
	}
	return(0);
}

int allocate_data_arrays(simulation_data& data_in, \
  const simulation_parameters& params_in) noexcept {
	const long max_simulation_count = \
	  (params_in.perform_start_number_sweep == 1) ? \
	  params_in.simulation_trials*params_in.sweep_count : \
	  params_in.simulation_trials;
	data_in.extinction_times = new long[max_simulation_count];
	fill(data_in.extinction_times, \
		  data_in.extinction_times + max_simulation_count, -1);
	if (params_in.keep_trajectory_samples == 1) {
		data_in.single_group_samples = \
		  (static_cast<long>(floor( \
		    static_cast<double>(params_in.max_timesteps)/ \
		    static_cast<double>(params_in.sample_interval) \
		)) + 1)*species_types*params_in.simulation_trials;
		const long max_samples = (params_in.perform_start_number_sweep == 1) ? \
		  params_in.sweep_count*data_in.single_group_samples : \
		  data_in.single_group_samples;
		data_in.samples = new long[max_samples];
		fill(data_in.samples, data_in.samples + max_samples, -1);
	}
	if (params_in.perform_start_number_sweep == 1) {
  		data_in.sweep_average_extinction_times = new long[max_simulation_count];
  		fill(data_in.sweep_average_extinction_times, \
		  data_in.sweep_average_extinction_times + max_simulation_count, -1);
	}
	if (params_in.keep_cycles_to_extinction == 1) {
		data_in.max_cycles_list = new long[max_simulation_count];
		fill(data_in.max_cycles_list, \
		  data_in.max_cycles_list + max_simulation_count, -1);
	}
	if (params_in.keep_extinction_time_histogram == 1) {
		data_in.extinction_time_histogram = \
		  new long[params_in.extinction_time_histogram_length];
		fill(data_in.extinction_time_histogram, \
		  data_in.extinction_time_histogram + \
		  params_in.extinction_time_histogram_length, 0);
	}
	return(0);
}


int write_data(const simulation_data& data, \
  const simulation_parameters& params) noexcept {
	ostringstream ss;
	ss.str(string());
	ss.clear();
	ss << "_mu" << params.mu << "_sigma" << params.sigma << "_lambda" << \
	  params.lambda << "_preystart" << params.prey_start_number << \
	  "_predatorstart" << params.predator_start_number;
	if (params.perform_start_number_sweep == 1) {
		const string sweep_species = (params.species_to_sweep == 0) ? \
		  "prey" : "predator";
		ss << "sweep" << sweep_species << "_schange" << params.sweep_change;
	}
	const string file_suffix = ss.str() + ".txt";
	ss.str(string());
	ss.clear();
	ofstream fp_out;
	int output_success = 0;
	if (params.perform_start_number_sweep == 1) {
		const long sweep_start_value = (params.species_to_sweep == 0) ? \
		  params.prey_start_number : params.predator_start_number;
		output_success += open_output_file("extinction_times" + file_suffix, \
		  fp_out, "extinction_times");
		fp_out << "sweep_population_startval,extinction_time" << endl;
		for (long i=0; i<params.sweep_count; ++i) {
			long offset = i*params.simulation_trials;
			long sweep_value = sweep_start_value + (i*params.sweep_change);
			for (long j=0; j<params.simulation_trials; ++j) {
				if (data.extinction_times[offset + j] > neg_one_long) {
					fp_out << sweep_value << "," << \
					  data.extinction_times[offset + j] << '\n';
				}
			}
		}
		fp_out.flush();
		output_success += close_output_file(fp_out, "extinction_times");
		if (params.keep_trajectory_samples == 1) {
			output_success += open_output_file( \
			  "population_trajectory" + file_suffix, fp_out, \
			  "population_trajectory");
			fp_out << "sweep_population_startval,prey,predators" <<  endl;
			for (long i=0; i<params.sweep_count; ++i) {
				long offset = i*data.single_group_samples;
				long sweep_value = sweep_start_value + (i*params.sweep_change);
				for (long j=0; j<data.single_group_samples; j+=2) {
					fp_out << sweep_value << "," << \
					  data.samples[offset + j] << "," << \
					  data.samples[offset + j + 1] << '\n';
				}
			}
			fp_out.flush();
			output_success += close_output_file(fp_out, \
			  "population_trajectory");
		}
		if (params.keep_cycles_to_extinction == 1) {
			output_success += open_output_file( \
			  "cycles_to_extinction" + file_suffix, fp_out, \
			  "cycles_to_extinction");
			fp_out << "sweep_population_startval,cycles_to_extinction" << endl;
			for (long i=0; i<params.sweep_count; ++i) {
				long offset = i*params.simulation_trials;
				long sweep_value = sweep_start_value + (i*params.sweep_change);
				for (long j=0; j<params.simulation_trials; ++j) {
					if (data.max_cycles_list[offset + j] > neg_one_long) {
						fp_out << sweep_value << "," << \
						  data.max_cycles_list[offset + j] << '\n';
					}
				}
			}
			fp_out.flush();
			output_success += close_output_file(fp_out, \
			  "cycles_to_extinction");
		}
	}
	else {
		output_success += open_output_file("extinction_times" + file_suffix, \
		  fp_out, "extinction_times");
		for (long j=0; j<params.simulation_trials; ++j) {
  			if (data.extinction_times[j] > 0) {
  				fp_out << data.extinction_times[j] << '\n';
  			}
  		}
		fp_out.flush();
		output_success += close_output_file(fp_out, "extinction_times");
		if (params.keep_trajectory_samples == 1) {
			output_success += open_output_file( \
			  "population_trajectory" + file_suffix, fp_out, \
			  "population_trajectory");
			fp_out << "prey,predators" <<  endl;
			for (long j=0; j<data.single_group_samples; j+=2) {
				fp_out << data.samples[j] << "," << \
				  data.samples[j + 1] << '\n';
			}
			fp_out.flush();
			output_success += close_output_file(fp_out, \
			  "population_trajectory");
		}
		if (params.keep_cycles_to_extinction == 1) {
			output_success += open_output_file( \
			  "cycles_to_extinction" + file_suffix, fp_out, \
			  "cycles_to_extinction");
			fp_out << "cycles_to_extinction" << endl;
			for (long j=0; j<params.simulation_trials; ++j) {
				if (data.max_cycles_list[j] > neg_one_long) {
					fp_out << data.max_cycles_list[j] << '\n';
				}
			}
			fp_out.flush();
			output_success += close_output_file(fp_out, \
			  "cycles_to_extinction");
		}
	}
	if (params.keep_extinction_time_histogram == 1) {
		output_success += open_output_file( \
		  "extinction_time_histogram" + file_suffix, fp_out, \
		  "extinction_time_histogram");
		for (long i=0; i<params.extinction_time_histogram_length; ++i) {
			fp_out << i*params.extinction_time_histogram_bin_size << " " \
			  << data.extinction_time_histogram[i] << '\n';
		}
		fp_out.flush();
		output_success += close_output_file(fp_out, \
		  "extinction_time_histogram");
	}
	if (output_success) {
		cout << "Error in writing data output." << endl;
		return(1);
	}
	return(0);
}

//The following block of functions assist with the simulation operation.
inline long direct_binomial(const double& success_prob, \
  const long& trials) noexcept {
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

inline void advance_species(long*const __restrict__ species_counts_in, \
  const double& mu_in, const double& sigma_in,
  const double& lambda_in) noexcept {
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

inline int quadrant(const long*const __restrict__ species_counts_in, \
  const long& equil_prey_in, const long& equil_predator_in) noexcept {
	//This function determines what "quadrant" the
	//predator prey system is in.
	int quadrant_out = -1;
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

inline long simulation_run_base(const simulation_parameters& params, \
  long*const __restrict__ species_counts) noexcept {
	const long timestep_limit = params.max_timesteps;
	long last_timestep = 0;
	for (long i=0; i<timestep_limit; ++i) {
		if (i%50 == 0) {
			cout << i << endl;
		}
		advance_species(species_counts, params.mu, params.sigma, \
		  params.lambda);
		for (int j=0; j<species_types; ++j) {
			if (species_counts[j] < 1) {
				species_counts[j] = 0;
				last_timestep = i + 1;
			}
		}
		if (last_timestep) {
			break;
		}
	}
	return(last_timestep);
}

void simulate_base(const simulation_parameters& params, \
  simulation_data& data) noexcept {
	const long sweep_limit = (params.perform_start_number_sweep == 1) ? \
	  params.sweep_count : 1;
	const long trial_limit = params.simulation_trials;
	for (long sweep=0; sweep<sweep_limit; ++sweep) {
		long population_start[2] = {
			params.prey_start_number,
			params.predator_start_number
		};
		if (params.perform_start_number_sweep == 1) {
			population_start[params.species_to_sweep] += \
			  sweep*params.sweep_change;
		}
		for (long trial=0; trial<trial_limit; ++trial) {
			cout << "trial = " << trial << endl;
			long species_counts[2] = {
				population_start[0],
				population_start[1]
			};
			long offset = sweep*trial_limit;
			long last_timestep = simulation_run_base(params, species_counts);
			if (last_timestep) {
				data.extinction_times[offset + trial] = last_timestep;
			}
		}
	}
	return;
}


int main(int argc, char** argv) {
	if(argc < 2){
		cout << "No input file specified." << endl;
		return(1);
	}

	//The following block reads in the input file parameters, and
	//checks that they make sense.
	const string inputfile = argv[1];
	ifstream fp_in;
	simulation_parameters params_input = {};
	const int input_success = read_input_parameters(inputfile, fp_in, \
	  params_input);
	if (input_success) {
		cout << "Reading input parameters failed." << endl;
		return(1);
	}
	if (fp_in.is_open()) {
		cout << "Input file failed to close." <<  endl;
		return(1);
	}

	const simulation_parameters params = {
		params_input.mu,
		params_input.sigma,
		params_input.lambda,
		params_input.equil_prey,
		params_input.equil_predator,
		params_input.prey_start_number,
		params_input.predator_start_number,
		params_input.max_timesteps,
		params_input.simulation_trials,
		params_input.keep_trajectory_samples,
		params_input.sample_interval,
		params_input.perform_start_number_sweep,
		params_input.species_to_sweep,
		params_input.sweep_count,
		params_input.sweep_change,
		params_input.keep_cycles_to_extinction,
		params_input.keep_extinction_time_histogram,
		params_input.extinction_time_histogram_length,
		params_input.extinction_time_histogram_bin_size
	};

	simulation_data data;
	allocate_data_arrays(data, params);

	//This initializes the random number generator.
	autoinit_randgen();

	simulate_base(params, data);

	write_data(data, params);


	data.~simulation_data();


/*
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

			else if (abs(previous_quadrant - quadrant) == two_int) {
				cout << "quadrant jump, previous = " << previous_quadrant << \
				  ",  current = " << quadrant <<  ",  timestep = " << \
				  i << endl;
			}

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
				cout << "cycle count = " << cycle_count << endl;
				max_cycles_list[trial] = cycle_count;
				if (keep_all_extinction_times == 1) {
					extinction_times[trial] = last_timestep;
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
	}
*/

	return(0);
}
