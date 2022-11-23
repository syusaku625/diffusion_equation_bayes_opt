#include"one_diffusion.hpp"

using namespace std;

void onedimensinal_diffusion::input_parameter()
{
    if(material == "F"){
        string str,base_label,label,disp;
        base_label = "/Fluid";
        label = base_label + "/numOfelement";
        if ( !tp.getInspectedValue(label,numOfelement)){
          cout << label << " is not set" << endl;
          exit(0);
        }
        label = base_label + "/dx";
        if ( !tp.getInspectedValue(label,element_length)){
          cout << label << " is not set" << endl;
          exit(0);
        }
        label = base_label + "/dt";
        if ( !tp.getInspectedValue(label,dt)){
          cout << label << " is not set" << endl;
          exit(0);
        }
        label = base_label + "/diffusion_coefficient";
        if ( !tp.getInspectedValue(label,diffusion_coefficient)){
          cout << label << " is not set" << endl;
          exit(0);
        }
        label = base_label + "/outputDir";
        if ( !tp.getInspectedValue(label,outputDir)){
          cout << label << " is not set" << endl;
          exit(0);
        }
        label = base_label + "/boundary_node";
        if ( !tp.getInspectedValue(label,boundary_node)){
          cout << label << " is not set" << endl;
          exit(0);
        }
        label = base_label + "/boundary_value";
        if ( !tp.getInspectedValue(label,boundary_value)){
          cout << label << " is not set" << endl;
          exit(0);
        }
    }

    if(material == "S"){
        string str,base_label,label,disp;
        base_label = "/Solid";
        label = base_label + "/numOfelement";
        if ( !tp.getInspectedValue(label,numOfelement)){
          cout << label << " is not set" << endl;
          exit(0);
        }
        label = base_label + "/dx";
        if ( !tp.getInspectedValue(label,element_length)){
          cout << label << " is not set" << endl;
          exit(0);
        }
        label = base_label + "/dt";
        if ( !tp.getInspectedValue(label,dt)){
          cout << label << " is not set" << endl;
          exit(0);
        }
        label = base_label + "/diffusion_coefficient";
        if ( !tp.getInspectedValue(label,diffusion_coefficient)){
          cout << label << " is not set" << endl;
          exit(0);
        }
        label = base_label + "/outputDir";
        if ( !tp.getInspectedValue(label,outputDir)){
          cout << label << " is not set" << endl;
          exit(0);
        }
        label = base_label + "/boundary_node";
        if ( !tp.getInspectedValue(label,boundary_node)){
          cout << label << " is not set" << endl;
          exit(0);
        }
        label = base_label + "/boundary_value";
        if ( !tp.getInspectedValue(label,boundary_value)){
          cout << label << " is not set" << endl;
          exit(0);
        }
    }
}

void onedimensinal_diffusion::initialize()
{
    element.resize(numOfelement);
    c.resize(numOfelement+1);
    phi.resize(numOfelement+1);
    x.resize(numOfelement+1);
    K.resize(numOfelement+1);
    mass.resize(numOfelement+1);
    mass_inv.resize(numOfelement+1);
    for(int i=0; i<K.size(); i++){
        K[i].resize(numOfelement+1);
        for(int j=0; j<K[i].size(); j++){
            K[i][j]=0.0;
        }
    }
    for(int i=0; i<mass.size(); i++){
        mass[i].resize(numOfelement+1);
        for(int j=0; j<mass[i].size(); j++){
            mass[i][j]=0.0;
        }
    }
    for(int i=0; i<mass_inv.size(); i++){
        mass_inv[i].resize(numOfelement+1);
        for(int j=0; j<mass_inv[i].size(); j++){
            mass_inv[i][j]=0.0;
        }
    }
    for(int i=0; i<element.size(); i++){
        element[i].resize(2);
        element[i][0]=i;
        element[i][1]=i+1;
    }

    for(int i=0; i<x.size(); i++){
        x[i]=element_length*i;
    }
}

void onedimensinal_diffusion::calc_mass_matrix()
{
    for(int i=0; i<element.size(); i++){
        mass[element[i][0]][element[i][0]]+=1.0/3.0*element_length;
        mass[element[i][0]][element[i][1]]+=1.0/6.0*element_length;
        mass[element[i][1]][element[i][0]]+=1.0/6.0*element_length;
        mass[element[i][1]][element[i][1]]+=1.0/3.0*element_length;
    }
    for(int i=0; i<mass.size(); i++){
        double sum = 0.0;
        for(int j=0; j<mass[i].size(); j++){
            sum += mass[i][j];
            mass[i][j]=0.0;
        }
        mass[i][i]=sum;
    }
    for(int i =0; i<mass_inv.size(); i++){
        mass_inv[i][i] = 1.0/mass[i][i];
    }
}

void onedimensinal_diffusion::calc_K_matrix()
{
    for(int i=0; i<element.size(); i++){
        K[element[i][0]][element[i][0]] += diffusion_coefficient*1.0/element_length*element_phi[i];
        K[element[i][0]][element[i][1]] += -diffusion_coefficient*1.0/element_length*element_phi[i];
        K[element[i][1]][element[i][0]] += -diffusion_coefficient*1.0/element_length*element_phi[i];
        K[element[i][1]][element[i][1]] += diffusion_coefficient*1.0/element_length*element_phi[i];
    }
}

void onedimensinal_diffusion::set_initial_boundary_node()
{
    b_n.push_back(boundary_node);
    i_b_n.push_back(boundary_value);
    c[b_n[0]]=i_b_n[0];
}

void onedimensinal_diffusion::boundary_setting(std::vector<double> boundary)
{
    c[b_n[0]]=i_b_n[0];
}

void onedimensinal_diffusion::time_step(std::vector<double> boundary, int time)
{
    if(material=="F") boundary_setting(boundary);
    vector<double> R(c.size(), 0.0);
    vector<double> Dc(c.size(),0.0);
    vector<double> DcR(c.size(), 0.0);
    vector<double> MDcR(c.size(),0.0);
    for(int i=0; i<K.size(); i++){
        for(int j=0; j<c.size(); j++){
            Dc[i] += K[i][j]*c[j];
        }
    }
    for(int i=0; i<mass.size(); i++){
        for(int j=0; j<boundary.size(); j++){
            R[i]+=mass[i][j]*boundary[j];
        }
    }
    for(int i=0; i<Dc.size(); i++){
        DcR[i] = Dc[i]-R[i];
    }
    for(int i=0; i<mass.size(); i++){
        for(int j=0; j<DcR.size(); j++){
            MDcR[i] += mass_inv[i][j]*DcR[j];
        }
    }
    for(int i=0; i<c.size(); i++){
        c[i]=c[i]-dt*MDcR[i];
    }
}

void onedimensinal_diffusion::dump(int step)
{
    mkdir(outputDir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);

    string filename = outputDir+"/"+to_string(step) + ".dat";
    ofstream ofs(filename);
    for(int i=0; i<c.size(); i++){
        ofs << c[i] << endl;
    }
    ofs.close();
}

double onedimensinal_diffusion::access_c(int i)
{
    return c[i];
}

void onedimensinal_diffusion::setting_phi(double p)
{
    for(int i=0; i<phi.size(); i++) phi[i] = p;
}

double onedimensinal_diffusion::return_sum()
{
    double sum =0.0;
    for(int i=0; i<c.size(); i++){
        sum += c[i];
    }
    return sum;
}
