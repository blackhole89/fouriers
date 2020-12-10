#include <stdio.h>
#include <armadillo>
#include <vector>

using namespace arma;

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

int main(int argc, char* argv[])
{
    constexpr int N=4;
    fmat F(1<<N,Poly<N>::max_long());
    for(int i=0;i<Poly<N>::max_long();++i) {
        Poly<N> p;
        if(p.from_long(i)) {
            for(int j=0;j<(1<<N);++j) {
                bool v = p.eval(j);
                F(j,i)=v?-1:1;
            }
        } else {
            for(int j=0;j<(1<<N);++j) {
                F(j,i)=0;
            }
        }
        printf("%s\n",p.to_str().c_str());
    }
    F.print("F:");

    fvec And(1<<N,fill::ones);
    And((1<<N)-1)=-1;

    And.print("And:");

    fvec Coeffs = solve(F,And);
    Coeffs.print("Coeffs:");
    return 0;
}
