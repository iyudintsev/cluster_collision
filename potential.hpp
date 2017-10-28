#ifndef potential_hpp
#define potential_hpp

#include "Eigen/Eigen/Dense"
#include <iostream>

using namespace std;
using namespace Eigen;


class LennardJonesModel{
	double e_energy, e_force, g, g2;
	double r_min, r_min2, r_max2;
public:
	LennardJonesModel(){};
	LennardJonesModel(double e_val, double g_val);
	
	double calcEnergy(Vector3d r);
	Vector3d calcForce(Vector3d r);
};


#endif