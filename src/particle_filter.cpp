/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;
static default_random_engine gen;
void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	num_particles = 101;
	
	//Generate Gussian noise
	default_random_engine generator;
	normal_distribution<double> x_ran(0, std[0]);
	normal_distribution<double> y_ran(0, std[1]);
	normal_distribution<double> theta_ran(0, std[2]);
	//Generate particles
	for (int i=0; i<num_particles; i++){
		Particle p;
		p.id = i;
		p.x = x;
		p.y = y;
		p.theta = theta;
		p.weight = 1.0;
		//Add Guassian noise
		double random_x = x_ran(generator);
		double random_y = y_ran(generator);
		double random_theta = theta_ran(generator);
		p.x += random_x;
		p.y += random_y;
		p.theta += random_theta;

		particles.push_back(p);
	}
	is_initialized = true;
	
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	
	//Random generator
	double x_f;
	double y_f;
	double theta_f;
	default_random_engine generator;
	normal_distribution<double> x_ran(0, std_pos[0]);
	normal_distribution<double> y_ran(0, std_pos[1]);
	normal_distribution<double> theta_ran(0, std_pos[2]);
	//Loop over particles
	for (int i=0; i<num_particles; i++){
		//Predection steps
		if (fabs(yaw_rate)<0.00001){
			x_f = particles[i].x + velocity*delta_t*cos(particles[i].theta);
			y_f = particles[i].y + velocity*delta_t*sin(particles[i].theta);
			theta_f = 0.0;
		}
		else{
			x_f = particles[i].x + velocity/yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
			y_f = particles[i].y + velocity/yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
			theta_f = particles[i].theta + yaw_rate*delta_t;			
		}

		//Generate Guassian num for x, y, theta
		//Assign prediction with Guassian noise
		double random_x = x_ran(generator);
		particles[i].x = x_f + random_x;
		double random_y = y_ran(generator);
		particles[i].y = y_f + random_y;
		double random_theta = theta_ran(generator);
		particles[i].theta = theta_f + random_theta;		
	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	
	for (unsigned int i=0; i<observations.size(); i++){
		int map_id = -1;
		float min_dist = numeric_limits<float>::max();
		for (unsigned int j=0; j<predicted.size(); j++){
			float cal_dist = dist(observations[i].x, observations[i].y,
			predicted[j].x, predicted[j].y);
			if (cal_dist < min_dist){
				min_dist = cal_dist;
				map_id = predicted[j].id;
			}
		}
	observations[i].id = map_id;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html
	
	//Define landmark standard deviations;
	double sigma_x = std_landmark[0];
	double sigma_y = std_landmark[1];
	
	for (int i=0; i<num_particles; i++){
		//Define and assign particle coordinates in map coordinates system
		double rot_ang = particles[i].theta;
		double xp = particles[i].x;
		double yp = particles[i].y;
		//Created vector to store landmarks that are predicted to be in sensor range
		vector<LandmarkObs> predictions;
		//For each landmark, find its associated observaion
		for (unsigned int j=0; j<map_landmarks.landmark_list.size(); j++){
			double lm_x = map_landmarks.landmark_list[j].x_f;
			double lm_y = map_landmarks.landmark_list[j].y_f;
			int lm_id = map_landmarks.landmark_list[j].id_i;
			if (dist(lm_x, lm_y, xp, yp) < sensor_range){
				predictions.push_back(LandmarkObs{lm_id, lm_x, lm_y});
			}
		}
		//Create a copy of observations 
		vector<LandmarkObs> Obs_landmark;
		for (unsigned int k=0; k<observations.size(); k++){
			//Transform particles Obs_landmark coordiantes into map coordinates system
			double xc = observations[k].x;
			double yc = observations[k].y;
			double trans_x = xp + cos(rot_ang)*xc - sin(rot_ang)*yc;
			double trans_y = yp + sin(rot_ang)*xc + cos(rot_ang)*yc;
			Obs_landmark.push_back(LandmarkObs{observations[k].id, trans_x, trans_y});
		}
		//Associate observations with landmarks observed within certain range
		dataAssociation(predictions, Obs_landmark);
		//reinit weight
		particles[i].weight = 1.0;
		//Find each observation that associated with landmark and calculate the its map coordinates
		for (unsigned int j=0; j<Obs_landmark.size(); j++){
			//Define landmark coordinates in map coordinates system
			double mu_x;
			double mu_y;
			//Define Obs_landmark coordiantes in map coordiantes system
			double t_x = Obs_landmark[j].x;
			double t_y = Obs_landmark[j].y;			
			for (unsigned int k=0; k<predictions.size(); k++){
				if (Obs_landmark[j].id == predictions[k].id){
					//Get landmark coordinates (Map coordinates system)
					mu_x = predictions[k].x;
					mu_y = predictions[k].y;
				}
			}
			//Product of the observations weights is particle weight
			double power_term = -(pow((t_x-mu_x),2)/(2*pow(sigma_x,2)) + pow((t_y-mu_y),2)/(2*pow(sigma_y,2)));
			particles[i].weight *= 1/(2*M_PI*std_landmark[0]*std_landmark[1])*exp(power_term);
		}
	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
/* 	vector<Particle> new_particles;
	vector<double> weights;
	//Get current weight
	for (int i=0; i< particles.size(); i++){
		weights.push_back(particles[i].weight);
		cout<<"particles " << i << "weight is "<< particles[i].weight << endl;
	}
	//Generate random stargin index
	default_random_engine gen;
	uniform_int_distribution<int> uniintdist(0, num_particles-1);
	int index = uniintdist(gen);

	//Wheel spining
	//get maximum weight
	double max_weight = *max_element(weights.begin(), weights.end());
	//Generator for random value between 0.0 and max_weight
	uniform_real_distribution<double> unirealdist(0.0, max_weight);
	double beta = 0.0;
	for (int i=0; i<num_particles; i++){
		beta = beta + unirealdist(gen)*2.0;
		while (beta > weights[index]){
			beta -= weights[index];
			index = (index+1)%num_particles;
		}
		new_particles.push_back(particles[index]);
	}
	particles = new_particles; */
	  vector<Particle> new_particles;

  // get all of the current weights
  vector<double> weights;
  for (int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
  }

  // generate random starting index for resampling wheel
  uniform_int_distribution<int> uniintdist(0, num_particles-1);
  auto index = uniintdist(gen);

  // get max weight
  double max_weight = *max_element(weights.begin(), weights.end());

  // uniform random distribution [0.0, max_weight)
  uniform_real_distribution<double> unirealdist(0.0, max_weight);

  double beta = 0.0;

  // spin the resample wheel!
  for (int i = 0; i < num_particles; i++) {
    beta += unirealdist(gen) * 2.0;
    while (beta > weights[index]) {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    new_particles.push_back(particles[index]);
  }

  particles = new_particles;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
