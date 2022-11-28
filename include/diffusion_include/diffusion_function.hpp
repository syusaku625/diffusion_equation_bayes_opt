/*
-------------------------------------------------------------------------
   This file is part of BayesOpt, an efficient C++ library for 
   Bayesian optimization.

   Copyright (C) 2011-2015 Ruben Martinez-Cantin <rmcantin@unizar.es>
 
   BayesOpt is free software: you can redistribute it and/or modify it 
   under the terms of the GNU Affero General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   BayesOpt is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Affero General Public License for more details.

   You should have received a copy of the GNU Affero General Public License
   along with BayesOpt.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------
*/

#ifndef __TEST_FUNCTIONS_HPP__
#define __TEST_FUNCTIONS_HPP__

#define _USE_MATH_DEFINES
#include <H5Cpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include<string>
#include <boost/math/constants/constants.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include "bayesopt/bayesopt.hpp"
#include "specialtypes.hpp"
#include "two_diffusion.hpp"

class DiffusionModel: public bayesopt::ContinuousModel
{
public:
    DiffusionModel(bayesopt::Parameters par):
    ContinuousModel(6,par) {}

    void Init(std::string input_file){
        Fluid.phase_set("F");
        Solid.phase_set("S");
        Vessel.phase_set("V");
        std::cout << "input infomation" << std::endl;
        Fluid.input_info(input_file);
        Solid.input_info(input_file);
        Vessel.input_info(input_file);

        std::cout << "gauss_point_setting" << std::endl;
        Fluid.gauss_point_setting();
        Solid.gauss_point_setting();
        Vessel.gauss_point_setting();

        std::cout << "matrix initialize" << std::endl;
        Fluid.matrix_initialize();
        Solid.matrix_initialize();
        Vessel.matrix_initialize();

        std::cout << "calc matrix" << std::endl;
        Fluid.calc_matrix();
        Solid.calc_matrix();

        std::cout << "set boundary" << std::endl;
        Q_vc.resize(Vessel.numOfNode); Q_cv.resize(Vessel.numOfNode);
        Q_vi.resize(Vessel.numOfNode); Q_iv.resize(Vessel.numOfNode);
        Q_ci.resize(Vessel.numOfNode); Q_ic.resize(Vessel.numOfNode);
        Vessel.boundary_setting(0.0, Q_cv, Q_iv);

        C_sum.resize(Vessel.numOfElm);
        input_evaluate_phi(Vessel.numOfElm);
        iter_count=0;
    };

    double evaluateSample( const vectord& xin){
        
        for(int i=0; i<Vessel.time; i++){
            std::vector<double> element_C_vessel(Vessel.numOfElm), element_C_Fluid(Fluid.numOfElm), element_C_Solid(Solid.numOfElm);
            Vessel.transform_point_data_to_cell_data(element_C_vessel, Vessel.C);
            Fluid.transform_point_data_to_cell_data(element_C_Fluid, Fluid.C);
            Solid.transform_point_data_to_cell_data(element_C_Solid, Solid.C);
            for(int j=0; j<Vessel.numOfElm; j++){
                //vessel to csf
                if(element_C_vessel[j]>=element_C_Fluid[j]){
                  Q_vc[j] = xin(0)*sqrt(Vessel.phi[j])*(element_C_vessel[j]-element_C_Fluid[j]);
                }
                //csf to vessel
                else if(element_C_vessel[j]<element_C_Fluid[j]){
                  Q_vc[j] = xin(1)*sqrt(Vessel.phi[j])*(element_C_vessel[j]-element_C_Fluid[j]);
                }
                //vessel to isf
                if(element_C_vessel[j]>=element_C_Solid[j]){
                  Q_vi[j] = xin(2)*sqrt(Vessel.phi[j])*(element_C_vessel[j]-element_C_Solid[j]);
                }
                //isf to vessel
                else if(element_C_vessel[j]<element_C_Solid[j]){
                  Q_vi[j] = xin(3)*sqrt(Vessel.phi[j])*(element_C_vessel[j]-element_C_Solid[j]);
                }
                //csf to isf
                if(element_C_Fluid[j]>=element_C_Solid[j]){
                  Q_ci[j] = xin(4)*sqrt(Fluid.phi[j])*(element_C_Fluid[j]-element_C_Solid[j]);
                }
                //isf to csf
                else if(element_C_Fluid[j]<element_C_Solid[j]){
                  Q_ci[j] = xin(5)*sqrt(Fluid.phi[j])*(element_C_Fluid[j]-element_C_Solid[j]);
                }
            }
            Vessel.boundary_setting(Vessel.dt*i, Q_vc, Q_vi);
            Fluid.time_step(Q_vc, Q_ci, Vessel.dt*i);
            Solid.time_step(Q_vi, Q_ci, Vessel.dt*i);

            if(i%Vessel.output_interval==0){
                std::vector<double> Vessel_phiC(Vessel.numOfElm),CSF_phiC(Fluid.numOfElm),ISF_phiC(Solid.numOfElm);
                Vessel.transform_point_data_to_cell_data_phi(Vessel_phiC, Vessel.C);
                Fluid.transform_point_data_to_cell_data_phi(CSF_phiC, Fluid.C);
                Solid.transform_point_data_to_cell_data_phi(ISF_phiC, Solid.C);
                for(int j=0; j<Vessel.numOfElm; j++){
                  C_sum[j] = Vessel_phiC[j] + CSF_phiC[j] + ISF_phiC[j];
                }
            }
        }
        
        double evaluation=0.0;
        for(int i=0; i<Vessel.numOfElm; i++){
            evaluation += pow(evaluation_phi[i]-C_sum[i],2.0);
        }
        evaluation /= Vessel.numOfElm;
        std::cout << "r_vc:" << xin(0) << " " << "r_cv:" << xin(1) << " " << "r_vi:" << xin(2) << " " << "r_iv:" << xin(3) << " " << "r_ci:" << xin(4) << " " << "r_ic:" << xin(5) << " " << "J:" << evaluation << std::endl;
        hdf5_dump("r_parameter.h5", iter_count, xin(0), xin(1), xin(2), xin(3), xin(4), xin(5 ));

        Vessel.reset();
        Fluid.reset();
        Solid.reset();
        iter_count++;
        return evaluation;
    };

  bool checkReachability(const vectord &query)
  {return true;};

  inline double sqr( double x ){ return x*x; };

  void printOptimal()
  {
    vectord sv(6);
    sv <<= 0.20169, 0.150011, 0.476874, 0.275332, 0.311652, 0.6573;

    std::cout << "Solution: " << sv << "->" 
	      << evaluateSample(sv) << std::endl;
  }

  private:
    twodimensinal_diffusion Solid, Vessel, Fluid;
    std::vector<double> C_sum;
    std::vector<double> Q_vc;
    std::vector<double> Q_cv;
    std::vector<double> Q_vi;
    std::vector<double> Q_iv;
    std::vector<double> Q_ci;
    std::vector<double> Q_ic;
    std::vector<double> evaluation_phi;
    double iter_count=0;
    double evaluation;

    void input_evaluate_phi(int numOfElm)
    {
      evaluation_phi.resize(numOfElm);
      std::ifstream ifs("evaluation_C.dat");
      if(!ifs){
          std::cout << "can't open evaluation_C.dat" << std::endl;
      }
      std::string str;
      for(int i=0; i<Vessel.numOfElm; i++){
          getline(ifs,str);
          evaluation_phi[i] = stod(str);
      }
      ifs.close();
    }
    
    void hdf5_export_parameter_and_cost_function(H5::H5File &file, const std::string &dataName, const double i_data);
    void hdf5_dump(const std::string output_h5_name, const int ic, const double r_vc, const double r_cv, const double r_vi, const double r_iv, const double r_ci, const double r_ic);


};


#endif
