#include "model.hpp"
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <random>


void Model::initRefl(int L0){
	reflection.push_back(Vector3d::Zero());
	for(int i = 0; i < 3; i++){
		for (int j = -1; j < 2; j++){
			if (j == 0) continue;
			Vector3d refl(0,0,0);
			refl(i) = j * L0;
			reflection.push_back(refl);
		}
	}
}


Model::Model(double L0, int cell_num0, double h0, int num_threads0){
	L = L0;
	initRefl(L0);

	h = h0;
	h2_2 = h * h / 2;
	h_2 = .5 * h;

	N = 0;
	cell_num = cell_num0;
	cell_len = L / cell_num;

	double g = 1;
	lg = LennardJonesModel(1, g);
	min_dist = pow(2, 1./6) * g;

	data_file = "./view/data.txt";

	distances = new int[4];
	distances[0] = 0;
	distances[1] = 1;
	distances[2] = cell_num;
	distances[3] = cell_num * cell_num;

	num_threads = num_threads0;
	if (num_threads > 0){
		Nth = N / num_threads;
	}
}


Model::~Model(){
	delete[] distances;
}


void Model::createParticle(Vector3d r, Vector3d v, double m){
	r = periodicBoundCond(r);
	Particle p(r, v, m, N);
	particles.push_back(p);
	N += 1;
}


Vector3d Model::periodicBoundCond(Vector3d r){
	for (int i = 0; i < 3; i++){
		if (r(i) < 0){
			while(r(i) < 0){
				r(i) += L;
			}
		}
		else if (r[i] > L){
			while(r(i) > L){
				r(i) -= L;
			}	
		} 
		else continue;
	}
	return r;
}


void Model::updateNeighbors(){
	vector<int> indices;
	for (int i = 0; i < N; i++){
		int cell = findCell(particles[i].r);
		indices.push_back(cell);
	}

	for (int i = 0; i < N; i++){
		Particle & p = particles[i];
		p.clearNeighbors();
		int index = indices[i];
		for (int j = 0; j < N; j++){
			if (i == j){
				continue;
			}
			for (int k = 0; k < 4; k++){
				if (abs(indices[j] - index) == distances[k]){
					p.addNeighbor(&particles[j]);
					break;
				}
			}
		}
	}
}


int Model::findCell(Vector3d r){
	int index[3];
	r = periodicBoundCond(r);
	for (int i = 0; i < 3; i++){
		index[i] = int(r(i) / cell_len) + 1;
	}
	return index[0] + index[1] * cell_num + index[2] * cell_num * cell_num;
}


void Model::init(int n, const vector<Vector3d>& coordinates,
			     const vector<Vector3d>& velocity, 
			     const vector<double>& mass){
	for (int i = 0; i < n; i++){
		Vector3d r = coordinates[i];
		r = periodicBoundCond(r);
		createParticle(r, velocity[i], mass[i]);
	}
	updateNeighbors();

	ofstream dump_file;
	dump_file.open(data_file);
	dump_file.close();
}


void Model::updateCoordinate(Particle & p){
	p.r += p.v * h + p.f * h2_2;
}


void Model::calcForcesWithThreads(int ibeg, int iend){
	for(int i = ibeg; i < iend; i++){
		Particle & pi = particles[i];
		for(Particle * pj : pi.neighbors){
			Vector3d dr = pi.r - pj->r;
			for (Vector3d& refl: reflection){
				Vector3d force = lg.calcForce(dr+refl);
				pi.f += force;
			}
		}
	}
	cout << "done\n";
}

void Model::calcForcesWithoutThreads(){
	for(int i = 0; i < N; i++){
		Particle & pi = particles[i];
		for(Particle * pj : pi.neighbors){
			Vector3d dr = pi.r - pj->r;
			for (Vector3d& refl: reflection){
				Vector3d force = lg.calcForce(dr+refl);
				pi.f += force;
			}
		}
	}
}


void Model::calcForces(){
	for(Particle& p: particles){
		p.f.setZero();
	}
	if (num_threads > 0){
		thread workers[num_threads];
		for(int i=0; i < num_threads; i++){
			int ibeg = i*Nth;
			int iend = (i == num_threads - 1) ? N : (i+1)*Nth;
			workers[i] = thread(&Model::calcForcesWithThreads, this, ibeg, iend);
		}

		for(int i=0; i < num_threads; i++){
			workers[i].join();
		}
	}
	else{
		calcForcesWithoutThreads();
	}
}


void Model::updateSpeed(Particle & p, const Vector3d & f0){
	p.v += (p.f + f0) * h_2;
}


void Model::algStep(){
	vector<Vector3d> forces;
	for(Particle & p: particles){
		updateCoordinate(p);
		forces.push_back(p.f);
	}
	calcForces();
	for (int i = 0; i < N; i++){
		Particle & p = particles[i];
		updateSpeed(p, forces[i]);
	}
}


void Model::run(int number_of_steps, int neighbors_count, int dump_num){
	for (int num_step = 1; num_step <= number_of_steps; num_step++){
		if (num_step % neighbors_count == 0){
			updateNeighbors();
		}
		if (num_step % dump_num == 0){
			cout << "Step: " << num_step << endl;
			dumpCoordinates();
		}

		algStep();
		for (Particle& p: particles){
			p.r = periodicBoundCond(p.r);
		}
	}
}


void Model::dumpCoordinates(){
	ofstream dump_file;
	dump_file.open(data_file, fstream::app);
	for (Particle& p: particles){
		dump_file << p.r.transpose() << endl;
	}
	dump_file << "next" << endl;
	dump_file.close();
}
