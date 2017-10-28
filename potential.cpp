#include "potential.hpp"


LennardJonesModel::LennardJonesModel(double e_val, double g_val){ 
	g = g_val;
	g2 = g_val * g_val;
	e_energy = 4 * e_val;
	e_force = 24 * e_val / g2;
	r_min = 0.88 * g;
	r_min2 = r_min * r_min;
	double r_max = 2.5 * g;
	r_max2 = r_max * r_max;
};


double LennardJonesModel::calcEnergy(Vector3d r){
	double r2 = r.transpose() * r;
	double u = 0;
	
	if (r2 > r_max2){
		return 0;
	}

	if (r2 < r_min2){
		Vector3d dr = r;
		double r_norm = r.norm();
		r = (r_min / r_norm) * r;
		dr = dr - r;
		u -= dr.transpose() * calcForce(r);
		r_norm = r.norm();
		r2 = r_norm * r_norm;
	}

	double ratio2 = g2 / r2;
	double ratio6 = ratio2 * ratio2 * ratio2;
	double ratio12 = ratio6 * ratio6;
	u += e_energy * (ratio12 - ratio6);	
	return u;
};


Vector3d LennardJonesModel::calcForce(Vector3d r){
	double r2 = r.transpose() * r;

	if(r2 > r_max2){
		return Vector3d::Zero();
	}
	
	if(r2 < r_min2){
		double r_norm = r.norm();
		r = (r_min / r_norm) * r;
		r_norm = r.norm();
		r2 =  r_norm * r_norm;
	}
	double ratio2 = g2 / r2;
	double ratio4 = ratio2 * ratio2;
	double ratio8 = ratio4 * ratio4;
	double ratio14 = ratio8 * ratio4 * ratio2;
	return e_force * (2 * ratio14 - ratio8) * r;
};
