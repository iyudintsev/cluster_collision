#include "model.hpp"
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <ctime>
#include <random>
#define MASTER 0
#define FROM_MASTER 1
#define FROM_WORKER 2

MPI_Status status;


void Model::initRefl(int l_size){
	reflection.push_back(Vector3d::Zero());
	for(int i = 0; i < 3; i++){
		for (int j = -1; j < 2; j++){
			if (j == 0) continue;
			Vector3d refl(0,0,0);
			refl(i) = j * l_size;
			reflection.push_back(refl);
		}
	}
}


Model::Model(double l_size, int cell_number, double step){
	L = l_size;
	initRefl(l_size);

	h = step;
	h2_2 = h * h / 2;
	h_2 = .5 * h;

	N = 0;
	cell_num = cell_number;
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
}


Model::~Model(){
	delete[] distances;
	delete[] borders;
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
}


void Model::convertCoordinates(double * coor){
	for (int i=0; i<N; i++){
		Particle & pi = particles[i];
		for (int j=0; j<3; j++){
			pi.r(j) = coor[3*i + j];
		}
	}
}


void Model::convertForces(double * forces, int worker){
	int i0 = borders[worker-1];
	for (int i=i0; i<borders[worker]; i++){
		Particle & pi = particles[i];
		for (int j=0; j<3; j++){
			pi.f(j) = forces[3 * (i - i0)+j];
		}
	}
}


void Model::updateCoordinate(Particle & p){
	p.r += p.v * h + p.f * h2_2;
}


void Model::calcForces(int ibeg, int iend){
	for(Particle& p: particles){
		p.f.setZero();
	}
	for(int i = ibeg; i < iend; i++){
		Particle & pi = particles[i];
		for(int j = 0; j < N; j++){
			if (i == j) continue;
			Particle & pj = particles[j];
			Vector3d dr = pi.r - pj.r;
			for (Vector3d& refl: reflection){
				Vector3d force = lg.calcForce(dr+refl);
				pi.f += force;
			}
		}
	}
}


void Model::updateSpeed(Particle & p, const Vector3d & f0){
	p.v += (p.f + f0) * h_2;
}


void Model::run(int number_of_steps, int dump_num){
	time_t t0, t1;
	t0 = time(NULL);

	int taskid, numtasks;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
	MPI_Comm_size(MPI_COMM_WORLD,&numtasks);
	int numworkers = numtasks-1;

	int Nth = N / numworkers;
	borders = new int[numworkers + 1];
	borders[numworkers] = N;
	for (int i=0; i < numworkers; i++){
		borders[i] = i * Nth;
	}

	for (int num_step = 1; num_step <= number_of_steps; num_step++){
		vector<Vector3d> forces;
		double coor[3*N];

		if (taskid == MASTER){
			if (num_step % dump_num == 0){
				cout << "Step: " << num_step << endl;
			}

			for(Particle & p: particles){
				updateCoordinate(p);
				p.r = periodicBoundCond(p.r);
				forces.push_back(p.f);
			}

			for (int i=0; i<N; i++){
				Particle & pi = particles[i];
				for (int j=0; j<3; j++){
					coor[3*i+j] = pi.r(j);
				}
			}

			for (int worker=1; worker <= numworkers; worker++){
				MPI_Send(&coor, 3*N, MPI_DOUBLE, worker, FROM_MASTER, MPI_COMM_WORLD);
			}

			for (int worker=1; worker <= numworkers; worker++){
				int n = borders[worker]	- borders[worker-1];
				double lforces[3*n];
				MPI_Recv(&lforces, 3*n, MPI_DOUBLE, worker, FROM_WORKER, MPI_COMM_WORLD, &status);
				convertForces(lforces, worker);						
			}			
			
			for (int i = 0; i < N; i++){
				Particle & p = particles[i];
				updateSpeed(p, forces[i]);
			}
		}
		
		if (taskid != MASTER){
			MPI_Recv(&coor, 3*N, MPI_DOUBLE, MASTER, FROM_MASTER, MPI_COMM_WORLD, &status);
			convertCoordinates(coor);
			calcForces(borders[taskid-1], borders[taskid]);

			int n = borders[taskid]	- borders[taskid-1];
			double lforces[3*n];
			int i0 = borders[taskid-1];
			for (int i=i0; i<borders[taskid]; i++){
				Particle & pi = particles[i];
				for (int j=0; j<3; j++)
					lforces[3 * (i - i0) + j] = pi.f(j);
			}
			MPI_Send(&lforces, 3*n, MPI_DOUBLE, MASTER, FROM_WORKER, MPI_COMM_WORLD);
		}
	}

	if (taskid == MASTER){
		t1 = time(NULL);
		cout << "t=" << t1-t0 << endl;
	}

	MPI_Finalize();
}
