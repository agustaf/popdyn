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

//These are compile-time constant definitions for use in the simulation.
constexpr int species_types = 2;
constexpr int prey_int = 0;
constexpr int predator_int = 1;
constexpr int false_int = 0;
constexpr int true_int = 1;
constexpr long use_equil_vals_long = -1;
constexpr double zero_double = 0.0;
constexpr double one_double = 1.0;
constexpr int zero_int = 0;
constexpr int one_int = 1;
constexpr long zero_long = 0;
constexpr long one_long = 1;
constexpr long neg_one_long = -1;

typedef enum sweep_choice {
	no_sweep = 0,
	prey_start_number_sweep = 1,
	predator_start_number_sweep = 2,
	mu_sweep = 3,
	sigma_sweep = 4,
	lambda_sweep = 5
} sweep_choice;

sweep_choice sweep_parameter_from_string( \
  const string& sweep_parameter_string) noexcept {
	sweep_choice sweep_parameter = no_sweep;
	if (sweep_parameter_string == "prey_start_number") {
		sweep_parameter = prey_start_number_sweep;
	}
	else if (sweep_parameter_string == "predator_start_number") {
		sweep_parameter = predator_start_number_sweep;
	}
	else if (sweep_parameter_string == "mu") {
		sweep_parameter = mu_sweep;
	}
	else if (sweep_parameter_string == "sigma") {
		sweep_parameter = sigma_sweep;
	}
	else if (sweep_parameter_string == "lambda") {
		sweep_parameter = lambda_sweep;
	}
	else {
		cout << "Unknown sweep parameter string: " << \
		  sweep_parameter_string << endl;
	}
	return sweep_parameter;
}

//The following block of code contains the data structures for the simulation.
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
	int perform_parameter_sweep;
	sweep_choice schoice;
	long sweep_count;
	double sweep_change;
	int keep_cycles_to_extinction;
	int keep_extinction_time_histogram;
	long extinction_time_histogram_length;
	long extinction_time_histogram_bin_size;
} simulation_parameters;

typedef struct run_parameters {
	public:
	double mu;
	double sigma;
	double lambda;
	long max_timesteps;
	long sample_interval;
} run_parameters;

typedef struct simulation_data {
	//This struct holds the results data from the simulation.
	public:
	long* extinction_times;
	long* extinction_time_histogram;
	double* sweep_average_extinction_times;
	long* max_cycles_list;
	long* samples;
	long single_run_samples;
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
		single_run_samples = zero_long;
		single_group_samples = zero_long;
		return;
	}
} simulation_data;


//The following block of functions assist in reading the input
//file and checking that input values make sense.
int open_input_file(const string& filename, ifstream& fp_in, \
  const string& description) noexcept {
	//This function opens a file to read it,
	//and checks that it opened correctly.
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
	//This function closes a file which was opened for reading,
	//and checks that it closed correctly.
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
	//This function opens a file to write to it,
	//and check that it opened correctly.
	if (fp_out.is_open()) {
		cout << "Error, output file already open for " << description << endl;
		return(1);
	}
	fp_out.open(filename, ios::out);
	if (!fp_out.is_open()) {
		cout << "Error, output file failed to open for " << description << \
		  endl;
		return(1);
	}
	fp_out.clear();
	fp_out.seekp(0, ios::beg);
	return(0);
}

int close_output_file(ofstream& fp_out, const string& description) noexcept {
	//This function closes a file which was opened for writing,
	//and checks that it closed correctly.
	if (!fp_out.is_open()) {
		cout << "Error, output file already closed for " << description << \
		  endl;
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
	//This function checks that a double value is between 0 and 1,
	//as a probability value should be.
	if ((value <= zero_double) || (value > one_double)) {
		cout << "Input error, value out of range for double variable " << \
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

int sweep_parameter_string_check(const string& value, const string& name) \
  noexcept {
	//This function checks that a string is either "prey" or "predator".
	if ((value != "prey_start_number") && \
	  (value != "predator_start_number") && \
	  (value != "mu") && (value != "sigma") && \
	  (value != "lambda")) {
		cout << "Input error, invalid sweep parameter string choice for " << \
		  name << endl;
		return(1);
	}
	return(0);
}

int read_input_parameters(const string& filename, ifstream& fp_in, \
  simulation_parameters& params) noexcept {
	//This function reads the input parameter for the simulation,
	//and stores the input in a simulation_parameters data structure.
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
	if (params.prey_start_number == use_equil_vals_long) {
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
	if (params.predator_start_number == use_equil_vals_long) {
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
	input_success += find_descriptor(fp_in, "perform_parameter_sweep=");
	fp_in >> params.perform_parameter_sweep;
	input_success += binary_intcheck(params.perform_parameter_sweep, \
	  "perform_parameter_sweep");
	if (params.perform_parameter_sweep == true_int) {
		input_success += find_descriptor(fp_in, "sweep_parameter=");
		string parameter_to_sweep_input = "";
		fp_in >> parameter_to_sweep_input;
		input_success += sweep_parameter_string_check( \
		  parameter_to_sweep_input, "parameter_to_sweep");
		params.schoice = \
		  sweep_parameter_from_string(parameter_to_sweep_input);
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
	if (params.keep_trajectory_samples == true_int) {
		input_success += find_descriptor(fp_in, "sample_interval=");
		fp_in >> params.sample_interval;
		input_success += positive_longcheck(params.sample_interval, \
		  "sample_interval");
		if (params.sample_interval > params.max_timesteps) {
			cout << "Error, trajectory sample_interval greater than ";
			cout << "max_timesteps" << endl;
			input_success += 1;
			return(1);
		}
	}
	input_success += find_descriptor(fp_in, "keep_cycles_to_extinction=");
	fp_in >> params.keep_cycles_to_extinction;
	input_success += binary_intcheck(params.keep_cycles_to_extinction, \
	  "keep_cycles_to_extinction");
	input_success += find_descriptor(fp_in, "keep_extinction_time_histogram=");
	fp_in >> params.keep_extinction_time_histogram;
	input_success += binary_intcheck(params.keep_extinction_time_histogram, \
	  "keep_extinction_time_histogram");
	if (params.keep_extinction_time_histogram == true_int) {
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
		cout << "Error, required input parameters ";
		cout << "missing or out of bounds." << endl;
		return(1);
	}
	return(0);
}


//The following block of functions assist with allocating memory for
//results data, and writing results data to files.
int allocate_data_arrays(simulation_data& data_in, \
  const simulation_parameters& params_in) noexcept {
	//This function allocates data arrays for results data
	//within a simulation_data data structure.
	const long max_simulation_count = \
	  (params_in.perform_parameter_sweep == true_int) ? \
	  params_in.simulation_trials*params_in.sweep_count : \
	  params_in.simulation_trials;
	data_in.extinction_times = new long[max_simulation_count];
	fill(data_in.extinction_times, \
		  data_in.extinction_times + max_simulation_count, neg_one_long);
	if (params_in.keep_trajectory_samples == true_int) {
		const double sample_count_double = \
		  floor(static_cast<double>(params_in.max_timesteps)/ \
		  static_cast<double>(params_in.sample_interval));
		data_in.single_run_samples = \
		  species_types*(static_cast<long>(sample_count_double) + 1);
		data_in.single_group_samples = data_in.single_run_samples* \
		  params_in.simulation_trials;
		const long max_samples = \
		  (params_in.perform_parameter_sweep == true_int) ? \
		  params_in.sweep_count*data_in.single_group_samples : \
		  data_in.single_group_samples;
		data_in.samples = new long[max_samples];
		fill(data_in.samples, data_in.samples + max_samples, neg_one_long);
	}
	if (params_in.perform_parameter_sweep == true_int) {
		data_in.sweep_average_extinction_times = \
		  new double[params_in.sweep_count];
		fill(data_in.sweep_average_extinction_times, \
		  data_in.sweep_average_extinction_times + params_in.sweep_count, \
		  neg_one_long);
	}
	if (params_in.keep_cycles_to_extinction == true_int) {
		data_in.max_cycles_list = new long[max_simulation_count];
		fill(data_in.max_cycles_list, \
		  data_in.max_cycles_list + max_simulation_count, neg_one_long);
	}
	if (params_in.keep_extinction_time_histogram == true_int) {
		data_in.extinction_time_histogram = \
		  new long[params_in.extinction_time_histogram_length];
		fill(data_in.extinction_time_histogram, \
		  data_in.extinction_time_histogram + \
		  params_in.extinction_time_histogram_length, zero_long);
	}
	return(0);
}

int write_data(const simulation_data& data, \
  const simulation_parameters& params) noexcept {
	//This function writes results data from a simulation_data
	//data structure to files.
	ostringstream ss;
	ss.str(string());
	ss.clear();
	ss << "_mu" << params.mu << "_sigma" << params.sigma << "_lambda" << \
	  params.lambda << "_preystart" << params.prey_start_number << \
	  "_predatorstart" << params.predator_start_number;
	string sweep_string = "NA";
	if (params.perform_parameter_sweep == true_int) {
		switch (params.schoice) {
			case prey_start_number_sweep:
				sweep_string = "preystart";
				break;
			case predator_start_number_sweep:
				sweep_string = "predatorstart";
				break;
			case mu_sweep:
				sweep_string = "mu";
				break;
			case sigma_sweep:
				sweep_string = "sigma";
				break;
			case lambda_sweep:
				sweep_string = "lambda";
				break;
			default:
				cout << "write_data, unknown params.schoice value" << endl;
				break;
		}
		ss << "_sweep" << sweep_string << "_schange" << \
		  params.sweep_change;
	}
	const string file_suffix = ss.str() + ".txt";
	ss.str(string());
	ss.clear();
	ofstream fp_out;
	int output_success = 0;
	if (params.perform_parameter_sweep == true_int) {
		double sweep_start_value = 0.0;
		switch (params.schoice) {
			case prey_start_number_sweep:
				sweep_start_value = \
				  static_cast<double>(params.prey_start_number);
				break;
			case predator_start_number_sweep:
				sweep_start_value = \
				  static_cast<double>(params.predator_start_number);
				break;
			case mu_sweep:
				sweep_start_value = params.mu;
				break;
			case sigma_sweep:
				sweep_start_value = params.sigma;
				break;
			case lambda_sweep:
				sweep_start_value = params.lambda;
				break;
			default:
				cout << "write_data, unkown params.schoice value" << endl;
				break;
		}
		output_success += open_output_file("extinction_times" + file_suffix, \
		  fp_out, "extinction_times");
		fp_out << sweep_string << "_value,extinction_time" << endl;
		for (long i=0; i<params.sweep_count; ++i) {
			long offset = i*params.simulation_trials;
			double sweep_value = sweep_start_value + \
			  (static_cast<double>(i)*params.sweep_change);
			for (long j=0; j<params.simulation_trials; ++j) {
				if (data.extinction_times[offset + j] > neg_one_long) {
					fp_out << sweep_value << "," << \
					  data.extinction_times[offset + j] << '\n';
				}
			}
		}
		fp_out.flush();
		output_success += close_output_file(fp_out, "extinction_times");
		output_success += \
		  open_output_file("sweep_average_extinction_times" + file_suffix, \
		  fp_out, "sweep_average_extinction_times");
		fp_out << sweep_string << "_value,average_extinction_time" << endl;
		for (long i=0; i<params.sweep_count; ++i) {
			double sweep_value = sweep_start_value + \
			  (static_cast<double>(i)*params.sweep_change);
			fp_out << sweep_value << "," << \
			  data.sweep_average_extinction_times[i] << '\n';
		}
		fp_out.flush();
		output_success += \
		  close_output_file(fp_out, "sweep_average_extinction_times");
		if (params.keep_trajectory_samples == true_int) {
			output_success += open_output_file( \
			  "population_trajectory" + file_suffix, fp_out, \
			  "population_trajectory");
			fp_out << sweep_string << "_value,prey,predators" <<  endl;
			for (long i=0; i<params.sweep_count; ++i) {
				long offset = i*data.single_group_samples;
				double sweep_value = sweep_start_value + \
				  (static_cast<double>(i)*params.sweep_change);
				for (long j=0; j<data.single_group_samples; j+=2) {
					if (data.samples[offset + j] > neg_one_long) {
						fp_out << sweep_value << "," << \
						  data.samples[offset + j] << "," << \
						  data.samples[offset + j + 1] << '\n';
					}
				}
			}
			fp_out.flush();
			output_success += close_output_file(fp_out, \
			  "population_trajectory");
		}
		if (params.keep_cycles_to_extinction == true_int) {
			output_success += open_output_file( \
			  "cycles_to_extinction" + file_suffix, fp_out, \
			  "cycles_to_extinction");
			fp_out << sweep_string << "_value,cycles_to_extinction" << endl;
			for (long i=0; i<params.sweep_count; ++i) {
				long offset = i*params.simulation_trials;
				double sweep_value = sweep_start_value + \
				  (static_cast<double>(i)*params.sweep_change);
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
			if (data.extinction_times[j] > neg_one_long) {
				fp_out << data.extinction_times[j] << '\n';
			}
		}
		fp_out.flush();
		output_success += close_output_file(fp_out, "extinction_times");
		if (params.keep_trajectory_samples == true_int) {
			output_success += open_output_file( \
			  "population_trajectory" + file_suffix, fp_out, \
			  "population_trajectory");
			fp_out << "prey,predators" <<  endl;
			for (long j=0; j<data.single_group_samples; j+=2) {
				if (data.samples[j] > neg_one_long) {
					fp_out << data.samples[j] << "," << \
					  data.samples[j + 1] << '\n';
				}
			}
			fp_out.flush();
			output_success += close_output_file(fp_out, \
			  "population_trajectory");
		}
		if (params.keep_cycles_to_extinction == true_int) {
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
	if (params.keep_extinction_time_histogram == true_int) {
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
	long count = 0;
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
	const long mu_change = direct_binomial(mu_in, species_counts_in[prey_int]);
	const long sigma_change = \
	  direct_binomial(sigma_in, species_counts_in[predator_int]);
	const long species_product = \
	  species_counts_in[prey_int]*species_counts_in[predator_int];
	const long lambda_change = direct_binomial(lambda_in, species_product);
	species_counts_in[prey_int] += (mu_change - lambda_change);
	species_counts_in[predator_int] += (lambda_change - sigma_change);
	return;
}

inline int quadrant(const long*const __restrict__ species_counts_in, \
  const long& equil_prey_in, const long& equil_predator_in) noexcept {
	//This function determines what "quadrant" the
	//predator prey system is in.
	int quadrant_out = -1;
	if (species_counts_in[predator_int] > equil_predator_in) {
		if (species_counts_in[prey_int] > equil_prey_in) {
			quadrant_out = 2;
		}
		else {
			quadrant_out = 1;
		}
	}
	else{
		if (species_counts_in[prey_int] > equil_prey_in) {
			quadrant_out = 3;
		}
		else{
			quadrant_out = 0;
		}
	}
	return(quadrant_out);
}

inline long simulation_run_base(const run_parameters& params, \
  long*const __restrict__ species_counts) noexcept {
	const long timestep_limit = params.max_timesteps;
	long last_timestep = neg_one_long;
	for (long i=0; i<timestep_limit; ++i) {
		/*if (i%50 == 0) {
			cout << i << endl;
		}*/
		advance_species(species_counts, params.mu, params.sigma, \
		  params.lambda);
		if (species_counts[prey_int] < 1) {
			species_counts[prey_int] = 0;
			last_timestep = i + 1;
		}
		if (species_counts[predator_int] < 1) {
			species_counts[predator_int] = 0;
			last_timestep = i + 1;
		}
		if (last_timestep > neg_one_long) {
			break;
		}
	}
	return(last_timestep);
}

inline long simulation_run_keep_trajectory( \
  const run_parameters& params, simulation_data& data, \
  long*const __restrict__ species_counts, const long& run_index) noexcept {
	const long timestep_limit = params.max_timesteps;
	const long sample_interval_limit = params.sample_interval;
	const long offset = run_index*data.single_run_samples;
	long last_timestep = neg_one_long;
	long trigger_index = 0;
	long local_sample_count = 0;
	for (long i=0; i<timestep_limit; ++i) {
		advance_species(species_counts, params.mu, params.sigma, \
		  params.lambda);
		if (species_counts[prey_int] < 1) {
			species_counts[prey_int] = 0;
			last_timestep = i + 1;
		}
		if (species_counts[predator_int] < 1) {
			species_counts[predator_int] = 0;
			last_timestep = i + 1;
		}
		if (last_timestep > neg_one_long) {
			data.samples[offset + local_sample_count] = \
			  species_counts[prey_int];
			data.samples[offset + local_sample_count + 1] = \
			  species_counts[predator_int];
			break;
		}
		++trigger_index;
		if (trigger_index == sample_interval_limit) {
			data.samples[offset + local_sample_count] = \
			  species_counts[prey_int];
			data.samples[offset + local_sample_count + 1] = \
			  species_counts[predator_int];
			local_sample_count += 2;
			trigger_index = 0;
		}
	}
	return(last_timestep);
}

void simulate_base(const simulation_parameters& params, \
  simulation_data& data) noexcept {
	const long sweep_limit = \
	  (params.perform_parameter_sweep == true_int) ? \
	  params.sweep_count : 1;
	const long trial_limit = params.simulation_trials;
	for (long sweep=0; sweep<sweep_limit; ++sweep) {
		if (params.perform_parameter_sweep == true_int) {
			cout << "sweep = " << sweep << endl;
		}
		double average_extinction_time = 0.0;
		long sweep_extinction_count = 0;
		long population_start[2] = {
			params.prey_start_number,
			params.predator_start_number
		};
		run_parameters run_params = {
			params.mu,
			params.sigma,
			params.lambda,
			params.max_timesteps,
			params.sample_interval
		};
		if (params.perform_parameter_sweep == true_int) {
			switch (params.schoice) {
				case prey_start_number_sweep:
					population_start[prey_int] += \
					  sweep*static_cast<long>(params.sweep_change);
					break;
				case predator_start_number_sweep:
					population_start[predator_int] += \
					  sweep*static_cast<long>(params.sweep_change);
					break;
				case mu_sweep:
					run_params.mu += sweep*params.sweep_change;
					break;
				case sigma_sweep:
					run_params.sigma += sweep*params.sweep_change;
					break;
				case lambda_sweep:
					run_params.lambda += sweep*params.sweep_change;
					break;
				default:
					cout << "simulate_base, unknown params.schoice" << endl;
					break;
			}
		}
		for (long trial=0; trial<trial_limit; ++trial) {
			cout << "trial = " << trial << endl;
			long species_counts[2] = {
				population_start[prey_int],
				population_start[predator_int]
			};
			long offset = sweep*trial_limit;
			long last_timestep = \
			  simulation_run_base(run_params, species_counts);
			if (last_timestep > neg_one_long) {
				data.extinction_times[offset + trial] = last_timestep;
				average_extinction_time += static_cast<double>(last_timestep);
				++sweep_extinction_count;
			}
			else {
				cout << "Warning, no system extinction." << endl;
			}
		}
		if (params.perform_parameter_sweep == true_int) {
			average_extinction_time /= \
			  static_cast<double>(sweep_extinction_count);
			data.sweep_average_extinction_times[sweep] = \
			  average_extinction_time;
		}
	}
	return;
}

void simulate_keep_trajectory(const simulation_parameters& params, \
  simulation_data& data) noexcept {
	const long sweep_limit = \
	  (params.perform_parameter_sweep == true_int) ? \
	  params.sweep_count : 1;
	const long trial_limit = params.simulation_trials;
	for (long sweep=0; sweep<sweep_limit; ++sweep) {
		if (params.perform_parameter_sweep == true_int) {
			cout << "sweep = " << sweep << endl;
		}
		double average_extinction_time = 0.0;
		long sweep_extinction_count = 0;
		long population_start[2] = {
			params.prey_start_number,
			params.predator_start_number
		};
		run_parameters run_params = {
			params.mu,
			params.sigma,
			params.lambda,
			params.max_timesteps,
			params.sample_interval
		};
		if (params.perform_parameter_sweep == true_int) {
			switch (params.schoice) {
				case prey_start_number_sweep:
					population_start[prey_int] += \
					  sweep*static_cast<long>(params.sweep_change);
					break;
				case predator_start_number_sweep:
					population_start[predator_int] += \
					  sweep*static_cast<long>(params.sweep_change);
					break;
				case mu_sweep:
					run_params.mu += sweep*params.sweep_change;
					break;
				case sigma_sweep:
					run_params.sigma += sweep*params.sweep_change;
					break;
				case lambda_sweep:
					run_params.lambda += sweep*params.sweep_change;
					break;
				default:
					cout << "simulate_base, unknown params.schoice" << endl;
					break;
			}
		}
		for (long trial=0; trial<trial_limit; ++trial) {
			cout << "trial = " << trial << endl;
			long species_counts[2] = {
				population_start[prey_int],
				population_start[predator_int]
			};
			long run_index = sweep*trial_limit + trial;
			long last_timestep = simulation_run_keep_trajectory(run_params, \
			  data, species_counts, run_index);
			if (last_timestep > neg_one_long) {
				data.extinction_times[run_index] = last_timestep;
				average_extinction_time += static_cast<double>(last_timestep);
				++sweep_extinction_count;
			}
			else {
				cout << "Warning, no system extinction." << endl;
			}
		}
		if (params.perform_parameter_sweep == true_int) {
			average_extinction_time /= \
			  static_cast<double>(sweep_extinction_count);
			data.sweep_average_extinction_times[sweep] = \
			  average_extinction_time;
		}
	}
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

	//This casts to input parameters into an immutable data structure,
	//to be certain they can't change in the rest of program operation.
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
		params_input.perform_parameter_sweep,
		params_input.schoice,
		params_input.sweep_count,
		params_input.sweep_change,
		params_input.keep_cycles_to_extinction,
		params_input.keep_extinction_time_histogram,
		params_input.extinction_time_histogram_length,
		params_input.extinction_time_histogram_bin_size
	};

	//This declares and allocates the data structures for holding results.
	simulation_data data;
	allocate_data_arrays(data, params);

	//This initializes the random number generator.
	autoinit_randgen();

	//This performs the simulation.
	if (params.keep_trajectory_samples == true_int) {
		simulate_keep_trajectory(params, data);
	}
	else {
		simulate_base(params, data);
	}

	//This writes the results data to files.
	write_data(data, params);

	//This frees the memory holding the results data.
	data.~simulation_data();

	return(0);
}

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
