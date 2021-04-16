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
#include <iterator>
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
	double M1 = mass_proj;
	double M2 = mass_target;
	double M3 = 1.00727647;
	double M4 = 13.00335;
	double T1 = input_energy;

	double Q_value = 2.7217441;
	double scattered_angle = (scattering_angle)*(pi/180.0);

	double a = (M3+M4);
	double b = (-2)*sqrt(M1*M3*T1)*cos(scattered_angle);
	double c = (T1*(M1-M4)) - (M4*Q_value);

	double determinant = (b*b)-(4*a*c);
	double root_T3 = ((-b) + sqrt(determinant))/(2*a);

	return (root_T3*root_T3);
}

static void CS_integrator(double proj_energy, std::map<double, double> change_with_energy, std::map<double, double> change_with_angle, double incremented_CS[180]) {
	
	double differential_cs[180];
	double energy_const = change_with_energy[0.94329];

	for (int i =0; i<101; i++){
		differential_cs[i] = (change_with_energy[proj_energy]/energy_const) * change_with_angle[100.00];
	}

	for (int i = 100; i<=170; i=i+2){
		differential_cs[i] = (change_with_energy[proj_energy]/energy_const) * change_with_angle[(double)i];
	}
	for (int i=101; i<170; i=i+2){
		differential_cs[i] = abs(differential_cs[i+1]+differential_cs[i-1])/2.0;
	}
	for (int i=170; i< 180; i++){
		differential_cs[i] = (change_with_energy[proj_energy]/energy_const) * change_with_angle[170.000] ;
	}
	for (int y = 0; y< 180; y++){
		double cs = differential_cs[y];
		incremented_CS[y] = (cs*(pi/180)*(1E-31)*(2*pi)*(sin((y+1)*pi/180)));
	}
}

static double CDF_distributor (double central_angle, double cross_section_integrated[180], double atoms_per_slice, double probability) {
	//central angle is put in in degrees and CM frame.

	//Find the cumulative density function
	double cdf[2*angular_range];
	cdf[0] = cross_section_integrated[(int)central_angle - angular_range]*atoms_per_slice;
	for (int i=1;i<2*angular_range; i++){
		cdf[i] = (cross_section_integrated[(int)central_angle - angular_range]*atoms_per_slice) + cdf[i-1];
	}
	cdf;

	//Normalize cdf
	for (int i =0; i<2*angular_range; i++){
		cdf[i] = cdf[i]/cdf[(2*angular_range - 1 )];
	}
	cdf;

	//Find the angle corresponding to the sample
	int counter = 0;
	while( probability >= cdf[counter]) {
		counter = counter++;
	}
	counter = counter + central_angle - angular_range;

	//Returns CM frame angle
	return counter;
}

static void cs_interpolator(double incremented_CS[180], std::map<double, std::vector<double>> &map, double proj_energy){

  auto lower = map.lower_bound(proj_energy) == map.begin() ? map.begin() : 
  --( map.lower_bound(proj_energy));
  auto upper = map.upper_bound(proj_energy);

  for (int i= 0; i<180;i++){
	 incremented_CS[i] = ( lower->second[i] + (upper->second[i] - lower->second[i]) * double(proj_energy-lower->first)/fabs(upper->first-lower->first) );
  }

}
 
static double energy_loss_calc (double proj_energy){
	std::map<double, double> energy_loss;
	
	energy_loss[4.3] = 2.69E-3;

	energy_loss[4.0] = 2.83E-3;
	
	energy_loss[3.75] = 2.966E-3;

	energy_loss[3.5] = 3.118E-3;

	energy_loss[3.25] = 3.289E-3;

	energy_loss[3.0]= 3.39E-3;

	energy_loss[2.94] = 3.54E-3;

	energy_loss[2.75] = 3.707E-3;

	energy_loss[2.5] = 3.966E-3;

	energy_loss[2.25] = 4.273E-3;

	energy_loss[2.0] = 4.63E-3;

	energy_loss[1.8] = 4.913E-3;

	energy_loss[1.7] = 5.094E-3;

	energy_loss[1.6] = 5.33E-3;

	energy_loss[1.5] = 5.513E-3;

	energy_loss[1.4] = 5.757E-3;

	energy_loss[1.3] = 6.029E-3;

	energy_loss[1.2] = 6.337E-3;

	energy_loss[1.1] = 6.684E-3;

	energy_loss[1.0] = 7.082E-3;

	energy_loss[0.9] = 7.548E-3;

	energy_loss[0.8] = 8.085E-3;

	energy_loss[0.7] = 8.75E-3;

	energy_loss[0.6] = 9.541E-3;

	energy_loss[0.5] = 10.6E-3;

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

	//Inputs
	double mass_first = 12.0; //incident particle
	double mass_second = 2.01410178; //stationary target
	double initial_energy = 8.6; //energy of incident
	double charge_first = 6.0; //charge of incident
	double charge_second = 1.0; //charge of target

	//Target Characteristics
	int target_layer_number = 1000;
	double layer_thickness = 45E-9;; //use correct units
	double target_slice_density = 1E18;

	//Main parameters
	#define proj_number 2E10
	#define proj_success 300
	double bias_factor = 10;
	//Initiate in CM frame/radians
	double central_angle;


	//Calculation arrays
	double cross_section_loop[180];
	double increment_CS[180];
	double scattering_angle;

	//Final Output - to be made into spectrum
	double successful_scatters[proj_success][2];

	//Dealing with file inputs for fusion cs
	std::map<double, double> change_with_energy;
	std::map<double, double> change_with_angle;

	ifstream fusion_files;
	fusion_files.open("EnergyChange.txt", ios::in);
	for (int i=0;i<126;i++){
		double energy, cross_section;
		fusion_files >> energy >> cross_section;
		change_with_energy.insert(std::pair<double, double>(energy, cross_section));
	}
	fusion_files.close();

	ifstream fusion_files2;
	fusion_files2.open("AngleChange.txt", ios::in);
	for (int i=0;i<36;i++){
		double angle, cross_section;
		fusion_files2 >> angle >> cross_section;
		change_with_angle.insert(std::pair<double, double>(angle, cross_section));
	}
	fusion_files2.close();

	//Finding energy of scattered deuteron
	kinematics_calc(output_energy_first, output_energy_second, mass_first, mass_second, initial_energy);
	double proj_energy = output_energy_second[0]*1;//MeV

	//Creating storage map for cross-sections for several energies
	std::map<double, std::vector<double>> cross_section_storage;
	double cross_section_holder[180];

	for(auto iterator = change_with_energy.begin(); iterator != change_with_energy.end(); iterator++) {
		proj_energy = iterator->first;
		CS_integrator(proj_energy, change_with_energy, change_with_angle, cross_section_holder);

		std::vector<double> incremented_CS;
		for (int y = 0; y< 180; y++){
			incremented_CS.push_back(cross_section_holder[y]);
		}

		cross_section_storage.insert(std::pair<double, std::vector<double>>((proj_energy+1E-10), incremented_CS));
	}
	
	//Sample from uniform distribution
	auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
	std::mt19937 generator(seed);
	std::uniform_real_distribution<double> flat_distribution(0.0, (1.0));

	

	int total_count = 0;
	int success_count = 0;

	//Running deuteron through target layers
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

			double reaction_probability = reaction_cs*target_slice_density*bias_factor*1E6;

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
	energylab_file.open("Scattering angle and Scattering Energy - deuteron fusion - LAB.txt", ios::out);
	for (int i = 0; i<proj_success; i++){
		energylab_file <<successful_scatters[i][0] << "   " << successful_scatters[i][1] << "\n" ;
	}
	energylab_file << total_count;
	energylab_file.close();

}