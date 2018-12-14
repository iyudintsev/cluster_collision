#include <math.h>
#include "Eigen/Eigen/Dense"

using namespace Eigen;


class GenerateCluster{
    double size;
    int N;
    int h;
    int counter;
public:
    vector<Vector3d> coordinates;
    vector<Vector3d> velocity;
    vector<double> mass;
    GenerateCluster(int n_val, int h_val, double cell_size){
        N = 2 * n_val + 1;
        h = 2 * h_val + 1;
        size = .5 *pow(2, 0.5) * cell_size;
        counter = 0;
    };

    void run(Vector3d r0, Vector3d v, double m){
        for (int i = 0; i < N; i++){
            for (int j = 0; j < N; j++){
                for (int k = 0; k < h; k++){
                    if ((i + j + k) % 2 == 0){
                        Vector3d r = r0 + size * Vector3d(i, j, k);
                        coordinates.push_back(r);
                        velocity.push_back(v);
                        mass.push_back(m);
                        counter += 1;
                    }
                }
            }
        }
    }

    int get_size(){
        return counter;
    }
};
