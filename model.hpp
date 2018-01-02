#ifndef model_hpp
#define model_hpp
#include <vector>
#include <string>
#include <mpi.h>
#include "Eigen/Eigen/Dense"
#include "potential.hpp"

using namespace std;
using namespace Eigen;


class Particle{
	vector<Particle*> neighbors;
public:
	int id;
	double m;
	Vector3d r, r0, f, v;
	
	Particle(Vector3d r_val, Vector3d v_val, double mass, int ID){
		r = r_val;
		r0 = r_val;
		m = mass;
		v = v_val;
		f = Vector3d::Zero();
		id = ID;
	}

	void addNeighbor(Particle* p){
		neighbors.push_back(p);
	}

	void clearNeighbors(){
		neighbors.clear();
	}
};


class Model{
	double L, cell_len, 
		   h, h2_2, h_2;
	int N, cell_num;
	int * distances;
	LennardJonesModel lg;
	vector<Particle> particles;

	vector<Vector3d> reflection;
	string data_file;
	int * borders;
public:
	double min_dist;
	double Upot, Ukin;

	Model(double l_size, int cell_number, double step);
	~Model();

	void initRefl(int l_size);
	void createParticle(Vector3d r, Vector3d v, double m);
	Vector3d periodicBoundCond(Vector3d r);
	void updateNeighbors();
	int findCell(Vector3d r);

	void init(int n, const vector<Vector3d>& coordinates,
			  const vector<Vector3d>& velocity, 
			  const vector<double>& mass);
	
	void updateCoordinate(Particle & p);
	void convertCoordinates(double * coor);
	void convertForces(double * forces, int worker);
	void calcForces(int ibeg, int iend);
	void updateSpeed(Particle & p, const Vector3d & f0);
	void run(int number_of_steps, int dump_num);
	void dumpCoordinates();

	void calcEnergy();  // need to implement
};


#endif