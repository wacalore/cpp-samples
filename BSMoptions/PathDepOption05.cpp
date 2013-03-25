#include "PathDepOption05.h"
#include <cmath>
#include <iostream>

double N(double x) // the standard normal cumulative distribution function
{
   double gamma = 0.2316419;     double a1 = 0.319381530;
   double a2    =-0.356563782;   double a3 = 1.781477937;
   double a4    =-1.821255978;   double a5 = 1.330274429;
   double pi    = 4.0*atan(1.0); double k  = 1.0/(1.0+gamma*x);
   if (x>=0.0)
   {
      return 1.0-((((a5*k+a4)*k+a3)*k+a2)*k+a1)
                  *k*exp(-x*x/2.0)/sqrt(2.0*pi);
   }
   else return 1.0-N(-x);
}


void GetZ(Matrix& Z, SamplePath& S, Matrix C, Vector S0, double r, double dt)
{
    int m=S.size();
    int d=S[0].size();
    Matrix Cinv = Inverse(C);
    Vector sigma(d);
    for(int i=0;i<d;i++) sigma[i] = sqrt(C[i] ^ C[i]);

    Z[0] = Cinv*(log(S[0]/S0) - dt*(r + (-0.5)*sigma*sigma))/sqrt(dt);
    for(int j=1;j<m;j++)
    {
        Z[j] = Cinv*(log(S[j]/S[j-1]) - dt*(r + (-0.5)*sigma*sigma))/sqrt(dt);
    }
}


void Rescale(SamplePath& S, Matrix& CZ, Matrix C, Vector S0, double r, double dt)
{
   int m=S.size();
   int d=S0.size();
   Vector sigma(d);
   for(int i=0;i<d;i++) sigma[i] = sqrt(C[i] ^ C[i]);

   S[0] = S0*exp(dt*(r + (-0.5*sigma*sigma)) + sqrt(dt)*(C*CZ[0]));
   for (int j=1; j<m; j++)
   {
      S[j] = S[j-1]*exp(dt*(r + (-0.5*sigma*sigma) + sqrt(dt)*(C*CZ[j])));
   }
}


double PathDepOption::PriceByMC(BSModel Model, long N, double epsilon)
{
    double H=0.0; double Hsq=0.0;
    int d = Model.S0.size();
    SamplePath S(m);
    Matrix C = Model.C;

    Matrix CZ(d);
    for(int i=0;i<d;i++) CZ[i].resize(m);

    delta.resize(d); rho.resize(d); vega.resize(d); theta.resize(d); gamma.resize(d);
    Vector Hdelta(d), Hvega(d), Hrho(d), Htheta(d), Hgamma(d);


    for (int i=0; i<d; i++)
    {
        delta[i] = 0.0; Hdelta[i] = 0.0; Hvega[i] = 0.0; Hrho[i] = 0.0; Htheta[i] = 0.0;
    }

    for(long i=0; i<N; i++)
    {
        Model.GenerateSamplePath(T,m,S);

        GetZ(CZ, S, C, Model.S0, Model.r, T/m);
        H = (i*H + Payoff(S))/(i+1.0);
        Hsq = (i*Hsq + pow(Payoff(S),2.0))/(i+1.0);

        for (int j=0; j<d; j++)
        {
            Vector S0tmp = Model.S0;
            S0tmp[j] = (1.0+epsilon)*S0tmp[j];
            Matrix Cvega = C;
            Cvega[j] = (1.0 + epsilon)*C[j];


            Rescale(S, CZ, C, S0tmp, Model.r, T/m);
            Hdelta[j] = (i*Hdelta[j] + Payoff(S)) / (i+1.0);

            Rescale(S, CZ, C, Model.S0, Model.r*(1.0+epsilon), T/m);
            Hrho[j] = (i*Hrho[j] + Payoff(S)) / (i+1.0);

            Rescale(S, CZ, Cvega, Model.S0, Model.r, T/m);
            Hvega[j] = (i*Hvega[j] + Payoff(S)) / (i+1.0);

            Rescale(S, CZ, C, Model.S0, Model.r, T*(1.0+epsilon)/m);
            Htheta[j] = (i*Htheta[j] + Payoff(S)) / (i+1.0);

            S0tmp[j] = (1.0-epsilon)*Model.S0[j];
            Rescale(S, CZ, C, S0tmp, Model.r, T/m);
            Hgamma[j] = (i*Hgamma[j] + Payoff(S)) / (i+1.0);

            Rescale(S, CZ, C, Model.S0, Model.r, T/m);
        }
    }

    for(int i=0; i<d; i++)
    {
        delta[i] = exp(-Model.r*T) * ( Hdelta[i] - H ) / (epsilon*Model.S0[i]);
        rho[i] = (exp(-Model.r*T*(1.0+epsilon))*Hrho[i] - exp(-Model.r*T)*H)/(epsilon*Model.r);
        vega[i] = exp(-Model.r*T) * (Hvega[i] - H) / (epsilon*sqrt(C[i]^C[i]));
        theta[i] = -1.0*(exp(-Model.r*T*(1.0+epsilon))*Htheta[i] - exp(-Model.r*T)*H)/(epsilon*T);
        gamma[i] = exp(-Model.r*T)*(Hdelta[i]-2.0*H+Hgamma[i])/(Model.S0[i]*Model.S0[i]*epsilon*epsilon);

    }

    Price = exp(-Model.r*T)*H;
    PricingError = exp(-Model.r*T)*sqrt(Hsq-H*H)/sqrt(N-1.0);

    return Price;
}

double EurBasketCall::Payoff(SamplePath& S)
{
    int d = S[0].size();
    double basket_value = 0.0;

    for(int i=0; i<d; i++)
    {
        basket_value += S[0][i];
    }

    if (basket_value > K) { return basket_value - K;}
    return 0.0;
}

double EurBasketPut::Payoff(SamplePath& S)
{
    int d = S[0].size();
    double basket_value = 0.0;

    for(int i=0; i<d; i++)
    {
        basket_value += S[0][i];
    }

    if (basket_value < K) { return K - basket_value; }
    return 0.0;

}

double GroupedEurCall::Payoff(SamplePath &S)
{
    int d = S[0].size();
    double payoff = 0.0;

    for (int i=0; i<d; i++)
    {
        if (S[0][i] > K[i]) payoff += S[0][i] - K[i];
    }

    return payoff;
}

double GroupedEurPut::Payoff(SamplePath& S)
{
    int d = S[0].size();
    double payoff = 0.0;

    for (int i=0; i<d; i++)
    {
        if (S[0][i] < K[i]) payoff += K[i] - S[0][i];
    }

    return payoff;
}

double GroupedEurCall::PriceByBSFormula(BSModel model)
{
    int d = model.S0.size();
    double price = 0.0;

    for (int j=0; j<d; j++)
    {
        EurCall Option(T, K[j]);
        price += Option.PriceByBSFormula(model.S0[j], model.sigma[j], model.r);
    }

    return price;
}


double GroupedEurPut::PriceByBSFormula(BSModel model)
{
    int d = model.S0.size();
    double price = 0.0;

    for (int j=0; j<d; j++)
    {
        EurPut Option(T, K[j]);
        price += Option.PriceByBSFormula(model.S0[j], model.sigma[j], model.r);
    }

    return price;
}

double PathDepOption::PriceByVarRedMC (BSModel model,long N,double epsilon,PathDepOption &G)
{
   DifferenceOfOptions HminusG(T,m,this,&G);

   Price = HminusG.PriceByMC(model,N,epsilon) + G.PriceByBSFormula(model);

   PricingError = HminusG.PricingError;

   return Price;
}

double ArthmAsianCall::Payoff(SamplePath& S)
{
    double Ave=0.0;
    int d=S[0].size();
    Vector one(d);

    for (int i=0; i<d; i++) one[i]=1.0;
    for (int k=0; k<m; k++)
    {
        Ave=(k*Ave+(one^S[k]))/(k+1.0);
    }
    if (Ave<K) return 0.0;
    return Ave-K;
}

double ArthmAsianPut::Payoff(SamplePath& S)
{
    double Ave=0.0;
    int d=S[0].size();
    Vector one(d);

    for (int i=0; i<d; i++) one[i]=1.0;
    for (int k=0; k<m; k++)
    {
        Ave=(k*Ave+(one^S[k]))/(k+1.0);
    }
    if (Ave > K) return 0.0;
    return K - Ave;
}


double EurCall::PriceByBSFormula(double S0, double sigma, double r)
{
   double d_plus = (log(S0/K)+(r+0.5*pow(sigma,2.0))*T)
            /(sigma*sqrt(T));
   double d_minus =  d_plus-sigma*sqrt(T);
   return S0*N(d_plus)-K*exp(-r*T)*N(d_minus);
}

double EurPut::PriceByBSFormula(double S0, double sigma, double r)
{
   double d_plus = (log(S0/K)+(r+0.5*pow(sigma,2.0))*T)
            /(sigma*sqrt(T));
   double d_minus =  d_plus-sigma*sqrt(T);
   return K*exp(-r*T)*N(-d_minus) - S0*N(-d_plus);
}
