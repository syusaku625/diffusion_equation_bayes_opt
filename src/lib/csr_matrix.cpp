#include"two_diffusion.hpp"

using namespace std;

void twodimensinal_diffusion::calc_adjacent_nodes()
{
  int tmp;
  bool b1;

  for(int ic=0;ic<numOfNode;ic++){
    for(int i=0,n=ieb[ic].size();i<n;i++){
      for(int j=0;j<element[ieb[ic][i]].size();j++){
        tmp = element[ieb[ic][i]][j];
        b1 = true;
        for(int k=0,n=inb[ic].size();k<n;k++){
          if(inb[ic][k]==tmp){
            b1 = false;
            break;
          }
        }
        if(b1==false) continue;
        inb[ic].push_back(tmp);
      }
    }
    sort(inb[ic].begin(),inb[ic].end());
  }
}

void twodimensinal_diffusion::calc_adjacent_elements()
{
  int tmp;
  for(int ic=0;ic<numOfElm;ic++){
    for(int j=0;j<element[ic].size();j++){
      tmp = element[ic][j];
      ieb[tmp].push_back(ic);
    }
  }
}

void twodimensinal_diffusion::CSR_initialize(const vector<vector<int>> &inb,const int &numOfNode,const int &dim)
{
  ptr.resize(numOfNode*dim+1);

  CSR_ptr_initialize(inb,numOfNode,dim);
  index.resize(nnz);
  value.resize(nnz);

  CSR_index_initialize(inb,numOfNode,dim);
}

void twodimensinal_diffusion::CSR_ptr_initialize(const vector<vector<int>> &inb,const int &numOfNode,const int &dim)
{
  nnz = 0;
  for(int i=0;i<dim;i++){
    for(int ic=0;ic<numOfNode;ic++){
      ptr[ic+i*numOfNode] = nnz;
      nnz += inb[ic].size()*dim;
    }
  }
  ptr[dim*numOfNode] = nnz;
}

void twodimensinal_diffusion::CSR_index_initialize(const vector<vector<int>> &inb,const int &numOfNode,const int &dim)
{
  int tmp = 0;
  for(int dim2=0;dim2<dim;dim2++){
    for(int ic=0;ic<numOfNode;ic++){
      for(int k=0;k<dim;k++){
        for(int i=0;i<inb[ic].size();i++){
          index[tmp] = inb[ic][i]+k*numOfNode;
          tmp++;
        }
      }
    }
  }
}

void twodimensinal_diffusion::set_CSR_value1D(vector<vector<vector<double>>> &K,const vector<vector<int>> &element,const int &numOfNode,
                               const int &numOfElm,const vector<vector<int>> &inb)
{
  int tmp1,tmp2,tmp3;

  #pragma omp parallel for
  for(int ic=0;ic<nnz;ic++) value[ic] = 0e0;

  for(int ielm=0;ielm<numOfElm;ielm++){
    for(int p=0;p<element[ielm].size();p++){

      tmp1 = element[ielm][p];
      for(int q=0;q<element[ielm].size();q++){

        tmp2 = element[ielm][q];

        for(int i=ptr[tmp1];i<ptr[tmp1+1];i++){
          if(tmp2==index[i]){
            value[i]     += K[ielm][p][q];
            break;
          }
        }
      }
    }
  }
}