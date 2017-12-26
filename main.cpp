#include "model.hpp"
#include "generator.cpp"
#include <ctime>


void clusterCollision(){
	Model model(20, 5, 0.01, 4);
	
	GenerateCluster generator(4, 4, model.min_dist);
	generator.run(Vector3d(13, 9, 9), Vector3d(-0.5, 0, 0), 1);
	generator.run(Vector3d(3, 8, 8), Vector3d(0.5, 0, 0), 1);

	int size = generator.get_size();
	cout << "Number of particles: " << size << endl;
	model.init(size, generator.coordinates, generator.velocity, generator.mass);
	model.run(500, 20, 50);
}


int main(){
	double t0, t1;
	t0 = time(NULL);
	clusterCollision();
	t1 = time(NULL);
	cout << "dt = " << t1 - t0 << endl;
	return 0;
}
