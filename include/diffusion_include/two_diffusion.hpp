#ifndef _DIFFUSION_H_
#define _DIFFUSION_H_

#include <H5Cpp.h>
#include<vector>
#include <sys/stat.h>
#include "TextParser.h"
#include<string>
#include<sstream>
#include<fstream>
#include<iostream>
#include<omp.h>
#include<set>
#include "mkl.h"
#include <stdio.h>
#include <stdlib.h>

class twodimensinal_diffusion{
    private:
        
        int numOfBoundaryNode;
        std::string material_judge;
        int numofomp;
        int numOfUpdatePoint;
        std::string gauss_setting, outputDir;
        std::string node_file, element_file, boundary_file, phi_file;
        std::vector<double> boundary_value;
        
        std::set<int> boundary_node_judge;
        std::string output_h5_name;
    public:
        double dt;
        double diffusion_coefficient;
        int time, output_interval, numOfNode, numOfElm;
        double coupling_coefficient_vi; //vessel to isf
        double coupling_coefficient_iv; //isf to vessel
        double coupling_coefficient_vc; //vessel to csf
        double coupling_coefficient_cv; //csf to vessel
        double coupling_coefficient_ci; //csf to isf
        double coupling_coefficient_ic; //isf to csf
        std::vector<std::vector<double>> mass;
        std::vector<double> mass_centralization;
        std::vector<int> boundary_node;
        std::vector<int> update_point;
        std::vector<double> C;
        std::vector<double> phi, phi_v;
        std::vector<std::vector<double>> node;
        std::vector<std::vector<int>> element;
        std::vector<std::vector<double>> D;
        std::vector<std::vector<double>> gauss;
        TextParser tp;
        void phase_set(std::string mat){
            material_judge=mat;
        }
        void export_vtu(const std::string &file, std::string judge, std::vector<double> output_value);
        void input_info(std::string input_file);
        void input_phi();
        int CountNumbersOfTextLines(std::string &filePath);
        void read_geometry();
        void boundary_initialize();
        void calc_dxdr(int ic, std::vector<std::vector<double>> node, std::vector<std::vector<int>> element, std::vector<std::vector<double>> &dxdr, std::vector<std::vector<double>> dNdr);
        void calc_dNdx(std::vector<std::vector<double>> &dNdx, std::vector<std::vector<double>> dNdr, std::vector<std::vector<double>> drdx);
        void calc_inverse_matrix_2x2(std::vector<std::vector<double>> dxdr, std::vector<std::vector<double>> &drdx);
        void gauss_point_setting();
        void matrix_initialize();
        void calc_matrix();
        void boundary_setting(double time, std::vector<double> Q_cv, std::vector<double> Q_iv);
        void time_step(std::vector<double> Q1, std::vector<double> Q2, double time);

        void dump(int ic);
        void hdf5_dump(int ic);
        void exportHDF5_double_1D(H5::H5File &file, const std::string &dataName, std::vector<double> i_data, int i_dim);
        double access_c(int ic);
        void transform_point_data_to_cell_data(std::vector<double> &element_C, std::vector<double> C);
        void transform_point_data_to_cell_data_phi(std::vector<double> &phiC, std::vector<double> C);
        void MKL_matrix_product(const std::vector<std::vector<double>> A_r, const std::vector<double> B_r, std::vector<double> &C_r, int m, int k, int n);

};
#endif
