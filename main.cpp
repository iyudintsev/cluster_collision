#include "model.hpp"
#include "generator.cpp"
#include <ctime>
#include <fstream>
#include <cstdlib>


void clusterCollision(int num_threads){
	Model model(20, 5, 0.01, num_threads);
	
	GenerateCluster generator(4, 4, model.min_dist);
	generator.run(Vector3d(13, 9, 9), Vector3d(-0.5, 0, 0), 1);
	generator.run(Vector3d(3, 8, 8), Vector3d(0.5, 0, 0), 1);

	int size = generator.get_size();
	cout << "Number of particles: " << size << endl;
	cout << "Number of threads: " << num_threads << endl;
	model.init(size, generator.coordinates, generator.velocity, generator.mass);
	model.run(100, 20, 20);
}


void dumpResults(time_t t0, time_t t1, int num_threads){
	double dt = t1 - t0;
	cout << "dt = " << dt << endl;

	ofstream out;
	out.open("results.dat", std::ios::app);
    struct tm * now = localtime( & t0 );
    out << (now->tm_year + 1900) << '-' 
        << (now->tm_mon + 1) << '-'
        << now->tm_mday << ' ' 
        << now->tm_hour << ':'
        << now->tm_min  << ':'
        << now->tm_sec  << ' ';
    out << "threads=" << num_threads << ' ';
    out << "dt=" << dt << endl;
    out.close();
}


int main(int argc, char *argv[]){
	int num_threads = 0;
	if (argc > 1){
		num_threads = atoi(argv[1]);
	}
	time_t t0, t1;
	t0 = time(NULL);
	clusterCollision(num_threads);
	t1 = time(NULL);
	dumpResults(t0, t1, num_threads);
	return 0;
}
