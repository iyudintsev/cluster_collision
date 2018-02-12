#include "model.hpp"
#include "generator.cpp"


void clusterCollision(){
	Model model(20, 5, 0.01);
	
	GenerateCluster generator(4, 4, model.min_dist);
	generator.run(Vector3d(13, 9, 9), Vector3d(-2, 0, 0), 1);
	generator.run(Vector3d(3, 8, 8), Vector3d(2, 0, 0), 1);

	int size = generator.get_size();
	model.init(size, generator.coordinates, generator.velocity, generator.mass);
	model.run(50, 10);
}


int main(){
	clusterCollision();
	return 0;
}
