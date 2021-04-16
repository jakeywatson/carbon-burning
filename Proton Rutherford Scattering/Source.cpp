
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <random>
#include <chrono>
#include <stdlib.h>
#include <map>
#include <array>
#include <vector>
#include <algorithm>
using namespace std;

//Constants
#define pi  3.14159265359
#define k  8.991804694E9
#define e  1.60217662E-19
#define angular_range 5
/*
	MODULE 1: KINEMATICS CALCULATION

	This module calculates the kinetic energy of an incident beam particle after having undergone elastic Rutherford scattering by a static target particle.

	INPUTS:
	mass_first      -  double - mass of the incident particle, amu
	mass_second     -  double - mass of the target particle, amu
	initial_energy  -  double - kinetic energy of the incident particle, MeV
	charge_first    -  double - proton number of the incident particle
	charge_second   -  double - proton number of the target particle

	OUTPUTS:
	output_energy_first - double array - array containing energy, in MeV, of incident particle after scattering, over 180 degrees
	output_energy_second - double array - array containing energy, in MeV of target particle after scattering, over 180 degrees
*/

static void kinematics_calc (double output_energy_first[180], double output_energy_second[180], double mass_first, double mass_second, double initial_energy, double scattering_angle = 0) {
	
	//Finds initial velocity of beam particle from input energy (non-relativistic) - sqrt(2mE)
	 double velocity_initial = sqrt(2*initial_energy/mass_first);

	//Determinant for the quadratic in final velocity.
	double determinant = (4*(mass_first*mass_first)*(velocity_initial*velocity_initial)*(cos(scattering_angle)*cos(scattering_angle))) - (4*(mass_first+mass_second)*(mass_first-mass_second)*(velocity_initial*velocity_initial));
	
	//Calculating maximum scattering angles for different mass cases (m1>m2, m2>m1)
	double projectile_max_angle_deflection_projectilegreater = (180/pi)*asin(mass_second/mass_first);
	double projectile_max_angle_deflection_targetgreater = 180;

	//Holders for velocities after scattering for incident and target particles
	double v_first;
	double v_second;

	//Checks to see if situation is m1 > m2 or m m2 < m1, applies maximum scattering angles accordingly. Fills output energy arrays for both particles.
	if (mass_second > mass_first) {

		//m2>m1: mass of target greater than mass of projectile

		//Filling ejectile energies
		for (int i = 0; i<projectile_max_angle_deflection_targetgreater; i++) {
			v_first = (2*mass_first*velocity_initial*cos(scattering_angle) + sqrt(determinant) )/(2*(mass_first+mass_second));
			output_energy_first[i] = 0.5*mass_first*v_first*v_first;

			//Incrementing angle
			scattering_angle = scattering_angle + (pi/180);
		}

		//Resetting scattering angle to zero for 2nd particle calculations
		scattering_angle = 0;

		//Filling recoil energies
		for (int i = 0; i<90; i++) {
			v_second = 2*(mass_first/(mass_second+mass_first))*velocity_initial*cos(scattering_angle);
			output_energy_second[i] = 0.5*mass_second*v_second*v_second;

			//Incrementing angle
			scattering_angle = scattering_angle + (pi/180);
		}
	} if (mass_first > mass_second) {
		
		//m1>m2: mass of projectile greater than mass of target

		//Filling incident particle energies after scattering
		for (int i = 0; i< projectile_max_angle_deflection_projectilegreater; i++) {
		
				v_first = (2*mass_first*velocity_initial*cos(scattering_angle) - sqrt(determinant) )/(2*(mass_first+mass_second));
				output_energy_first[i] = 0.5*mass_first*v_first*v_first;

				//Incrementing angle
				scattering_angle = scattering_angle + (pi/180);
			}

		//Resetting scattering angle to zero for 2nd particle calculations
		scattering_angle = 0;

		//Filling target particle energies after scattering
			for (int i = 0; i< 90; i++) {

				v_second = 2*(mass_first/(mass_second+mass_first))*velocity_initial*cos(scattering_angle);
				output_energy_second[i] = 0.5*mass_second*v_second*v_second;

				//Incrementing angle
				scattering_angle = scattering_angle +(pi/180);
			}
		}
}

static double final_energy_calc (double mass_proj, double mass_target, double input_energy, double scattering_angle){
	double final_velocity;
	scattering_angle = scattering_angle*(pi/180);

	//Finds initial velocity of beam particle from input energy (non-relativistic) - sqrt(2mE)
	 double velocity_initial = sqrt(2*input_energy/mass_proj);

	 double determinant = (4*(mass_proj*mass_proj)*(velocity_initial*velocity_initial)*(cos(scattering_angle)*cos(scattering_angle))) - (4*(mass_proj+mass_target)*(mass_proj-mass_target)*(velocity_initial*velocity_initial));

	 final_velocity = (2*mass_proj*velocity_initial*cos(scattering_angle) + sqrt(determinant) )/(2*(mass_proj+mass_target));

	 return (0.5*mass_proj*final_velocity*final_velocity);
}

/*
	MODULE 2: CROSS-SECTION CALCULATION

	This module calculates the Rutherford angular cross-section of an incident particle with a target particle.

	INPUTS:
	output_energy_second[i] -  double - energy of the particle scattered by incident 12C beam
	charge_first       -  double - proton number of the incident particle
	charge_second      -  double - proton number of the target particle

	OUTPUTS:
	cross_section      - double array - array containing angular cross section, for angles from 0-180, in increments of 1
*/

static void CS_calc (double cross_section[180], double particle_energy, double charge_first, double charge_second, double mass_first, double mass_second, double scattering_angle = 0) {
	
	double frame_correction = mass_first/(mass_first+mass_second);

	//Fills RF c-s array for increments of 1 degree, 1-180.
	for (int i = 0; i<180; i++) {
		
		cross_section[i] = pow( (k*charge_first*charge_second*e*e) / (4*(particle_energy*frame_correction)*e*pow(10, 6)), 2) * (1/(pow(sin(scattering_angle/2), 4)));
		
		//Incrementing angle
		scattering_angle = scattering_angle + (pi/180);
	}

	//Sets 0 degree c-s to constant, 1 degree c-s.
	cross_section[0] = cross_section[1];
}

/*
	MODULE 3: NUMERICAL INTEGRATOR

	This module calculates the total Rutherford cross-section of the incident particle with the target particle, through the method described at the start of the program.

	INPUTS:
	n_test_points - int - number of random points to be generated inside rectangle

	OUTPUTS:
	area - double - final value for the numerical integration, multiplied by 2 so as to account for both scattering directions. (0-2pi)
*/

static double CS_integrator (double cross_section[180], double n_test_points, double charge_first, double charge_second, double mass_first, double mass_second, double output_energy_second[180]) {
	double frame_correction = mass_first/(mass_first+mass_second);
	double area = 0;

	//Using high precision clock to give different seeds each time. Seeding mersenne twister RNG.
	auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
	std::mt19937 generator(seed);

	//Setting ranges for x, y generation
	std::uniform_real_distribution<double> x_uniform_distribution(0.0, (pi));
	std::uniform_real_distribution<double> y_uniform_distribution(0.0, cross_section[1]);

	double x_test, y_test, RF_predicted;
	int valid_points = 0;

		//Testing random points against predicted cs, discarding or saving
		for (int i = 0; i< n_test_points; i++) {
			 x_test = x_uniform_distribution(generator);
			 y_test = y_uniform_distribution(generator);

			 RF_predicted = pow( ((k*charge_first*charge_second*e*e) / (4*(output_energy_second[0]*frame_correction)*e*pow(10, 6))), 2) * (1/pow(sin(x_test/2), 4)) * sin(x_test);

			 if (RF_predicted >= y_test) {
				valid_points++;
			}
		}
		
		//Dividing number of saved points by total number generated, multiply by area of box of integration
		area = (((double)valid_points)/((double)n_test_points))*cross_section[1]*pi*pi*2;
		return (area);

	}

static void CS_integrator_increments (double cross_section_integrated[180], double cross_section[180], double n_test_points, double charge_first, double charge_second, double mass_first, double mass_second, double projectile_energy) {
	double frame_correction = mass_first/(mass_first+mass_second);
	double sum = 0;

	//Using high precision clock to give different seeds each time. Seeding mersenne twister RNG.
	auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
	std::mt19937 generator(seed);
	double area;

	for (int y=0; y< (180); y++) {
		
		double max = (y + 1.0)/(180);
		double min =(y + 0.0)/(180);
		//Setting ranges for x, y generation
		std::uniform_real_distribution<double> x_uniform_distribution(min, max);
		std::uniform_real_distribution<double> y_uniform_distribution(0.0, cross_section[y]);

		double x_test, y_test, RF_predicted;
		int valid_points = 0;

		//Testing random points against predicted cs, discarding or saving
		for (int i = 0; i< n_test_points; i++) {
			x_test = x_uniform_distribution(generator);
			 y_test = y_uniform_distribution(generator);

			 RF_predicted = pow( ((k*charge_first*charge_second*e*e) / (4*(projectile_energy*frame_correction)*e*pow(10, 6))), 2) * (1/pow(sin(x_test/2), 4)) * sin(x_test);

			 if (RF_predicted >= y_test) {
				valid_points++;
			}
		}

		//Dividing number of saved points by total number generated, multiply by area of box of integration
		area = (((double)valid_points)/((double)n_test_points))*cross_section[y]*(pi/180)*2*pi;
		//Adding to Array
		cross_section_integrated[y] = area;
	}

	/*
	Have issues with this section: if I miss out a factor of pi, answer to incremented cs agrees with non-incremented cs. Why is this?
	*/
}

static double CDF_distributor (double central_angle, double cross_section_integrated[180], double atoms_per_slice, double probability) {
	//central angle is put in in degrees and CM frame.

	//Find the cumulative density function
	double cdf[2*angular_range];
	cdf[0] = cross_section_integrated[(int)central_angle - angular_range]*atoms_per_slice;
	for (int i=1;i<2*angular_range; i++){
		cdf[i] = (cross_section_integrated[(int)central_angle - angular_range]*atoms_per_slice) + cdf[i-1];
	}

	//Normalize cdf
	for (int i =0; i<2*angular_range; i++){
		cdf[i] = cdf[i]/cdf[(2*angular_range - 1 )];
	}

	//Find the angle corresponding to the sample
	int counter = 0;
	while( probability >= cdf[counter]) {
		counter++;
	}
	counter = counter + central_angle - angular_range;

	//Returns CM frame angle
	return counter;
}

static void cs_interpolator(double incremented_CS[180], std::map<double, std::vector<double>> map, double proj_energy){

  auto lower = map.lower_bound(proj_energy) == map.begin() ? map.begin() : 
  --( map.lower_bound(proj_energy));
  auto upper = map.upper_bound(proj_energy);

  for (int i= 0; i<180;i++){
	 incremented_CS[i] = ( lower->second[i] + (upper->second[i] - lower->second[i]) * double(proj_energy-lower->first)/fabs(upper->first-lower->first) );
  }

}

static double energy_loss_calc (double proj_energy){
	std::map<double, double> energy_loss;

	energy_loss[3.0]= 2.01E-3;

	energy_loss[2.75] = 2.239E-3;

	energy_loss[2.5] = 2.403E-3;

	energy_loss[2.25] = 2.598E-3;

	energy_loss[2.0] = 2.32E-3;

	energy_loss[1.8] = 3.057E-3;

	energy_loss[1.7] = 3.186E-3;

	energy_loss[1.6] = 3.328-3;

	energy_loss[1.5] = 3.485E-3;

	energy_loss[1.4] = 3.660E-3;

	energy_loss[1.3] = 3.858E-3;

	energy_loss[1.2] = 4.083E-3;

	energy_loss[1.1] = 4.336E-3;

	energy_loss[1.0] = 4.593E-3;

	energy_loss[0.9] = 4.912E-3;

	energy_loss[0.8] = 5.293E-3;

	energy_loss[0.7] = 5.756-3;

	energy_loss[0.6] = 6.334E-3;

	energy_loss[0.5] = 7.082-3;

	auto lower = energy_loss.lower_bound(proj_energy) == energy_loss.begin() ? energy_loss.begin() : 
		--( energy_loss.lower_bound(proj_energy));
	auto upper = energy_loss.upper_bound(proj_energy);
	double energy_loss_value = ( lower->second + (upper->second - lower->second) * double(proj_energy-lower->first)/fabs(upper->first-lower->first) );

	return (proj_energy-energy_loss_value);
}

static double angle_of_interest(int i, double layer_thickness){
	double height_from_beam = 103.9305943E-3;
	double length_from_target = 138.9268929E-3;
	double angle = 180 - ((180/pi)*atan(height_from_beam/(length_from_target+(i*layer_thickness))));
	return angle;
}

int main() {
	
	//Outputs
	double output_energy_first[180];
	double output_energy_second[180];

	//Commented out the user input lines so as to streamline testing. Input values are hardcoded below, for a 6 MeV 12C incident on a target proton.
	/*
	cout << "Enter the mass of the first, incident particle (in amu): ";
	cin >> mass_first;
	cout << "Enter the mass of the second, target particle (in amu): ";
	cin >> mass_second;
	cout << "Enter the initial energy of the first, incident particle (in MeV): ";
	cin >> initial_energy;
	cout << "Enter the charge of the first particle: ";
	cin >> charge_first;
	cout << "Enter the charge of the second, target particle: ";
	cin >> charge_second;
	*/

	//Inputs
	double mass_first = 12.0; //incident particle
	double mass_second = 1.00727647; //stationary target
	double initial_energy = 8.6; //energy of incident
	double charge_first = 6.0; //charge of incident
	double charge_second = 1.0; //charge of target

	int target_layer_number = 1000;
	double layer_thickness = 45E-9;; //use correct units
	double target_slice_density = 1E18;

	#define proj_number 2000000
	#define proj_success 200
	double bias_factor = 5000;

	//Initiate in CM frame/radians
	double central_angle;

	//Calculation arrays
	double cross_section_loop[180];
	double increment_CS[180];
	double scattering_angle;

	//Final Output - to be made into spectrum
	double successful_scatters[proj_success][2];

	int success_counter = 0;
	//Sample from uniform distribution
	auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
	std::mt19937 generator(seed);
	std::uniform_real_distribution<double> flat_distribution(0.0, (1.0));

	//Finding energy of scattered proton
	kinematics_calc(output_energy_first, output_energy_second, mass_first, mass_second, initial_energy);
	double proj_energy = output_energy_second[0]*1;//MeV

	//Create cross-section array
	int div_number = (ceil((proj_energy-0.5)/0.1));
	std::map<double, std::vector<double>> cross_section_storage;

	for (int i = 0; i< div_number; i++) {
		CS_calc(cross_section_loop, proj_energy, charge_first, charge_second, mass_first, mass_second);

		CS_integrator_increments(increment_CS, cross_section_loop, 1000, charge_first, charge_second, mass_first, mass_second, proj_energy);

		std::vector<double> incremented_CS;
		for (int y = 0; y< 180; y++){
			incremented_CS.push_back(increment_CS[y]);
		}

		cross_section_storage.insert(std::pair<double, std::vector<double>>((proj_energy+1E-10), incremented_CS));
	
		proj_energy = proj_energy-0.1;
	}

	
	int total_count = 0;
	int success_count = 0;

	//Running proton through target layers
	for (int y=0; y< proj_number; y++){
		proj_energy = output_energy_second[0];
		total_count = total_count + 1;

		if (success_count >= proj_success){
			break;
		}

		for (int i = 0; i < target_layer_number; i++) {
			
			//cout << " \n Projectile: " << y << "  Layer: " << i << "\n";
			//cout << "Projectile Energy: " << proj_energy << "\n";
			double probability = flat_distribution(generator);
			double reaction_test = flat_distribution(generator);

			proj_energy = energy_loss_calc(proj_energy);
			central_angle = abs(pi - (2*(angle_of_interest(i, layer_thickness)*(pi/180))));

			//Integrate increments of cross-section
			cs_interpolator(increment_CS, cross_section_storage, proj_energy);

			//Explanation: cs_interpolator uses the current projectile energy, looks up and 
			//interpolates between the relevant cross-sections from the cross_section_storage map, and updates the current incremented_cs array.
			//Now can calculate reaction cs.

			double reaction_cs = 0;
			for (int z =0; z<(180); z++){
				(double)reaction_cs = (double)reaction_cs + increment_CS[z];
			}

			double reaction_probability = reaction_cs*target_slice_density*bias_factor;

			if (proj_energy > 0.5){

				if (reaction_test <= reaction_probability){
					
					double probability_detector_angles= 0;

					for (int i = 0; i< 2*angular_range; i++){
						int index = (central_angle*180/pi) - angular_range + i;
						probability_detector_angles = increment_CS[index] + probability_detector_angles;
					}
					probability_detector_angles = (probability_detector_angles/reaction_cs)*bias_factor;
					double scatter_test = flat_distribution(generator);

					if (scatter_test <= probability_detector_angles) {
						success_count++;

						//Find scattering angle CM FRAME - radians
						scattering_angle = CDF_distributor((central_angle*(180/pi)), cross_section_loop, target_slice_density, probability)*(pi/180);
			
						//Convert to LAB frame 
						scattering_angle = 0.5*(pi-scattering_angle);
				
						//Reconvert to degrees
						scattering_angle = 180 - scattering_angle*(180/pi);

						successful_scatters[success_count - 1][0] = scattering_angle; //LAB frame
						successful_scatters[success_count - 1][1] = final_energy_calc(mass_second, mass_first, proj_energy, scattering_angle); //In LAB frame
						//cout << "Successful " << success_count << "\n \n \n";
						break;

					} else {

					//cout << "Unsuccessful " << y+1 << "\n";
					//cout << "Success Count: " << success_count << "\n";
					break;

					}
				}

			} else {
				//If particle energy falls below 1 MeV, scrap and consider next particle
				break;
			}
		//Adding successful protons to array - angle and energy
		}
	}
	//Convolving data with Gaussian energy resolution of detector
	std::default_random_engine normal_generator;
	for (int i = 0; i<proj_success; i++){
		std:: normal_distribution<double> normal_distribution(successful_scatters[i][1], 0.0169864360);
		successful_scatters[i][1] = normal_distribution(normal_generator);
	}

	ofstream energylab_file;
	energylab_file.open("Scattering angle and Scattering Energy - proton 8.6 - LAB.txt", ios::out);
	for (int i = 0; i<proj_success; i++){
		energylab_file <<successful_scatters[i][0] << "   " << successful_scatters[i][1] << "\n" ;
	}
	energylab_file << total_count;
	energylab_file.close();
	////Writing output energies and cross section for both particles to file.
	//ofstream energylab_file;
	//energylab_file.open("Energy vs. Scattering angle -  Incident.txt", ios::out);
	//for (int i = 0; i<180; i++){
	//	energylab_file << output_energy_first[i] << "\n" ;
	//}
	//energylab_file.close();
	//ofstream energycm_file;
	//energycm_file.open("Energy vs. Scattering angle - Target.txt", ios::out);
	//for (int i = 0; i<180; i++){
	//	energycm_file << output_energy_second[i] << "\n" ;
	//}
	//energycm_file.close();
	//ofstream cs_file;
	//cs_file.open("Cross-Section vs. Scattering angle.txt", ios::out);
	//for (int i = 0; i<180; i++){
	//	cs_file << cross_section[i] << "\n" ;
	//}
	//cs_file.close();

	//Commented out the user input lines so as to streamline testing. Input values are hardcoded below, for a 6 MeV 12C incident on a target proton.
	/*
	cout << "Enter the mass of the first, incident particle (in amu): ";
	cin >> mass_first;
	cout << "Enter the mass of the second, target particle (in amu): ";
	cin >> mass_second;
	cout << "Enter the initial energy of the first, incident particle (in MeV): ";
	cin >> initial_energy;
	cout << "Enter the charge of the first particle: ";
	cin >> charge_first;
	cout << "Enter the charge of the second, target particle: ";
	cin >> charge_second;
	*/


	//exit
	return(0);
	}