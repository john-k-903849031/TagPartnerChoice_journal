#include <iostream>
#include <fstream>
#include "../../Empirical/include/emp/math/Random.hpp"
int main(){

	std::ofstream  file;
	file.open("means.csv");
	file << "mean,cutoff\n";
	emp::Random random(234);

	
	emp::vector<double> means = {-0.5, 0, 0.125, 0.1875, 0.25, 0.3125, 0.375, 0.4375, 0.5, 0.5625, 0.625, 1, 1.5, 2};
	
	for(int i = 0;i < 1000000; i++){
		for(auto mean : means){
			double cutoff = random.GetPoisson(mean*32);
			file << std::to_string(mean) << "," << std::to_string(cutoff) << "\n";
		}
	}
	file.close();
	return 0;
}
