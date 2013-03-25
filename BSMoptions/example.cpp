#include <iostream>
#include "PathDepOption05.h"


using namespace std;

int main()
{
   int d=3;
   Vector S0(d);
      S0[0]=40.0;
      S0[1]=60.0;
      S0[2]=100.0;
   double r=0.03;

   Matrix C(d);
   for(int i=0;i<d;i++) C[i].resize(d);


      C[0][0] =  0.1;  C[0][1] = -0.1;  C[0][2] = 0.0;
      C[1][0] = -0.1;  C[1][1] =  0.2;  C[1][2] = 0.0;
      C[2][0] =  0.0;  C[2][1] =  0.0;  C[2][2] = 0.3;

   BSModel Model(S0,r,C);

   double T=1.0/12.0, K=200.0;

   EurBasketPut Option(T,K);

   // initiating the control variate:
   double V=0.0;
   double epsilon = .003;
   for (int j=0;j<d;j++) V=V+Model.S0[j];
   Vector Kd=(K/V)*Model.S0;
   GroupedEurPut CVOption(T,Kd);

   long N=30000;
   Option.PriceByVarRedMC(Model,N, epsilon, CVOption);
   cout << "European Basket Put Price using Control Variates = "
        << Option.Price << endl
        << "Pricing Error = " << Option.PricingError << endl << endl;

   Option.PriceByMC(Model,N, epsilon);
   cout << "European Basket Put Price using direct MC        = "
        << Option.Price << endl
        << "Pricing Error = " << Option.PricingError << endl << endl;

   cout << "European Basket Put Greeks using Monte Carlo = " << endl;
   for(int i=0;i<d;i++) cout << "delta[" << i << "] = " << Option.delta[i]
                             << "\tgamma[" << i << "] = " << Option.gamma[i]
                             << "\tvega["  << i << "] = " << Option.vega[i] << endl;

   cout << "theta = " << Option.theta[0] << endl << "rho = " << Option.rho[0];
   return 0;
}

