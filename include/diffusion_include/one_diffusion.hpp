#ifndef _DIFFUSION_H_
#define _DIFFUSION_H_

#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include <sys/stat.h>
#include "TextParser.h"
#include<omp.h>

class onedimensinal_diffusion{
    private:
        std::string material;
        std::string outputDir;
        int boundary_node;
        double boundary_value;
        int numOfelement;
        double element_length;
        double dt;
        double diffusion_coefficient;
        std::vector<std::vector<int>> element;
        std::vector<double> c;
        std::vector<double> phi;
        std::vector<double> x;
        std::vector<std::vector<double>> mass;
        std::vector<std::vector<double>> mass_inv;
        std::vector<std::vector<double>> K;
        std::vector<int> b_n;
        std::vector<double> i_b_n;
    public:
        TextParser tp;
        std::vector<double> element_phi;
        void input_parameter();
        void initialize();
        void calc_mass_matrix();
        void calc_K_matrix();
        void set_initial_boundary_node();
        void time_step(std::vector<double> boundary, int time);
        void dump(int step);
        void setting_phi(double p);
        void boundary_setting(std::vector<double> boundary);
        double access_c(int i);
        double return_sum();

        onedimensinal_diffusion(std::string mat)
        {
            material = mat;
        }
};

#endif