#include<iostream>
#include<potential.hpp>
#include<fstream>
#include<iostream>

using namespace std;


class TestLennardJonesModel{
	double step, _step, error;
	LennardJonesModel lg;
public:
	TestLennardJonesModel(){
		lg = LennardJonesModel(1, 1);
		step = 1e-8;
		error = 1e-5;
		_step = 1. / step;
	};

	double scoreDif(const Vector3d & v, const Vector3d & vnum){
		return (v - vnum).norm() / v.norm();
	}

	void run(){
		bool success = true;
		for (int i=80; i<200; i++){
			Vector3d v(0.01 * i, 0, 0);
			Vector3d force = lg.calcForce(v);
			Vector3d num_force = numDerivative(v);
			double score = scoreDif(force, num_force);
			if (score > error){
				cout << "error: " << score << endl;
				cout << '\t' << force.transpose() << endl;
				cout << '\t' << num_force.transpose() << endl;
				success = false;
			};
		}
		if (success){
			cout << "Lennard-Jones Model: test was successful." << endl;
		}
	}

	void generatePlotData(char * file_name, bool force = false){
		ofstream data_file;
		data_file.open(file_name);

		for (int i=80; i<200; i++){
			Vector3d v(0.01 * i, 0, 0);
			double value;
			
			if(force){
				value = lg.calcForce(v)(0);
			} else value = lg.calcEnergy(v);

			data_file << v(0) << ' ' << value << endl;
		}

		data_file.close();
	}

	Vector3d numDerivative(const Vector3d & r){
		Vector3d num_force(0, 0, 0);
		double u = lg.calcEnergy(r);
		for (int i=0; i < 3; i++){
			Vector3d r_new = r;
			r_new(i) += step;
			double u_new = lg.calcEnergy(r_new);
			num_force(i) = (u_new - u) * _step;
		}
		return -num_force;
	}
};


int main(){
	TestLennardJonesModel test_lg;
	test_lg.run();
	// test_lg.generatePlotData("tests/data.txt");
	// test_lg.generatePlotData("tests/data.txt", true);
	return 0;
}