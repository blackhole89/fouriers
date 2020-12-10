#include <stdio.h>
#include <vector>
#include <numeric>

#include <coin/CoinBuild.hpp>
#include <coin/CoinModel.hpp>
#include <coin/ClpSimplex.hpp>
#include <coin/ClpMessage.hpp>
#include <coin/CoinMessageHandler.hpp>

#define itype int;

#define C2I(x,y) ((x)*N+(y))
#define AT(a,p) ((a)&(1<<(p)))

template<size_t N> struct Poly{
    std::vector<int> coeff;

    constexpr static int max_long() {
        return 1<<(N*(N+1)/2);
    }

    bool from_long(long coeffs) {
        bool is_deg2=false;;
        coeff.clear();
        coeff.resize(N*N);
        int c=0;
        for(int i=0;i<N;++i) {
            for(int j=i;j<N;++j) {
                if(AT(coeffs,c)) {
                    if(i!=j) is_deg2=true;
                    coeff[C2I(i,j)]=1;
                }
                ++c;
            }
        }
        return is_deg2;
    }

    bool eval(int input) {
        bool ret=false;
        for(int i=0;i<N;++i) {
            for(int j=0;j<N;++j) {
                if(coeff[C2I(i,j)] && AT(input,i) && AT(input,j)) ret=!ret;
            }
        }
        return ret;
    }

    std::string to_str() {
        std::string ret;
        for(int i=0;i<N;++i) {
            for(int j=0;j<N;++j) {
                if(coeff[C2I(i,j)]) {
                    if(ret.length()) ret+="+";
                    char buf[64];
                    if(i!=j) sprintf(buf,"x%dx%d",i,j);
                    else sprintf(buf,"x%d",i);
                    ret+=buf;
                }
            }
        }
        return ret;
    }
};

constexpr int N=4;

double evals(double* coeffs,std::vector<Poly<N>> &polies,int input)
{
    double ret=0.0;
    for(int i=0;i<Poly<N>::max_long();++i) {
        ret += (coeffs[2*i]-coeffs[2*i+1])*(polies[i].eval(input)?-1:1);
    }
    return ret;
}

int main(int argc, char* argv[])
{

    ClpSimplex model;
    model.setLogLevel(0);
    model.resize(0,2*Poly<N>::max_long());
    for(int i=0;i<Poly<N>::max_long()*2;++i) {
        model.setObjectiveCoefficient(i,1.0);
        model.setColumnLower(i,0.0);
        model.setColumnUpper(i,1.0);
    }

    std::vector<int> indices;
    std::vector<double> coeffs;
    indices.resize(2*Poly<N>::max_long());
    std::iota(indices.begin(),indices.end(),0);
    coeffs.resize(2*Poly<N>::max_long());

    std::vector<Poly<N> > polies;
    polies.resize(Poly<N>::max_long());
    for(int i=0;i<Poly<N>::max_long();++i) {
        polies[i].from_long(i);
        printf("%s\n",polies[i].to_str().c_str());
    }

    std::vector<int> And;
    And.resize(1<<N,1);
    And[(1<<N)-1]=-1;

    for(int j=0;j<(1<<N);++j) {
        for(int i=0;i<Poly<N>::max_long();++i) {
            bool v = polies[i].eval(j);
            coeffs[i*2]=v?-1:1;
            coeffs[i*2+1]=v?1:-1;
        }
        model.addRow(2*Poly<N>::max_long(), &indices[0], &coeffs[0], And[j], And[j]);
    }

    model.primal();
    
    double* sln = (double*)model.getColSolution();

    for(int i=0;i<Poly<N>::max_long();++i) {
        double v = sln[2*i]-sln[2*i+1];
        if(fabs(v)>0.00001) {
            printf("%.2f * %s\n",v,polies[i].to_str().c_str());
        }
    }
    printf("\n");
    for(int j=0;j<(1<<N);++j) {
        printf("@%d: %.2f. \n", j, evals(sln,polies,j));
    }

    return 0;
}
