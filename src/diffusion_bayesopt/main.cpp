#include <cmath>
#include <algorithm>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include "bayesopt/bayesopt.hpp"
#include<fstream>
#include "diffusion_function.hpp"

int main(int argc,char *argv[])
{
  //parameter settings
  bayesopt::Parameters par;
  par = initialize_parameters_to_default();
  par.n_iterations = 100;
  par.noise = 1e-2;
  par.random_seed = 0;
  par.verbose_level = 1;

  DiffusionModel diffusion6(par);
  diffusion6.Init(argv[1]);

  std::ofstream timelog;
  timelog.open("time_diffusion.log");
  std::clock_t curr_t;
  std::clock_t prev_t = clock();
  diffusion6.initializeOptimization();
  for (int i = 0; i < par.n_iterations; i++){      
    diffusion6.stepOptimization();
    curr_t = clock();
    timelog << i << "," << static_cast<double>(curr_t - prev_t) / CLOCKS_PER_SEC << std::endl;
    prev_t = curr_t;
  }

  timelog.close();

  vectord result = diffusion6.getFinalResult();
  std::cout << "Result: " << result << "->" << diffusion6.evaluateSample(result) << std::endl;

  return 0;
}