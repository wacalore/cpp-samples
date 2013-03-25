#ifndef PathDepOption05_h
#define PathDepOption05_h
#include "BSModel02.h"

class PathDepOption
{
    public:
        double T, Price, PricingError;
        int m;
        Vector delta, rho, vega, theta, gamma;

        double PriceByMC(BSModel Model, long N, double epsilon);
        virtual double PriceByBSFormula(BSModel model){return 0.0;}
        double PriceByVarRedMC(BSModel model,long N,double epsilon,PathDepOption &G);
        virtual double Payoff(SamplePath& S)=0;

};

class EurBasketCall: public PathDepOption
{
    public:
        double K;
        EurBasketCall(double T_, double K_) {T=T_; K=K_; m=1;}
        double Payoff(SamplePath &S);
};

class EurBasketPut: public PathDepOption
{
    public:
        double K;
        EurBasketPut(double T_, double K_) {T=T_; K=K_; m=1;}
        double Payoff(SamplePath &S);
};

class ArthmAsianCall: public PathDepOption
{
    public:
        double K;
        ArthmAsianCall(double T_, double K_, int m_)
        {T=T_; K=K_; m=m_;}
        double Payoff(SamplePath& S);
};

class ArthmAsianPut: public PathDepOption
{
    public:
        double K;
        ArthmAsianPut(double T_, double K_, int m_)
        {T=T_; K=K_; m=m_;}
        double Payoff(SamplePath& S);
};

class GroupedEurCall: public PathDepOption
{
    public:
        Vector K;
        GroupedEurCall(double T_, Vector K_) {T=T_; K=K_; m=1;}
        double PriceByBSFormula(BSModel model);
        double Payoff(SamplePath& S);
};

class GroupedEurPut: public PathDepOption
{
    public:
        Vector K;
        GroupedEurPut(double T_, Vector K_) {T=T_; K=K_; m=1;}
        double PriceByBSFormula(BSModel model);
        double Payoff(SamplePath& S);
};


class EurCall
{
    public:
      double K, T; // strike
      int m;
      EurCall(double T_, double K_)
         {T=T_; K=K_; m=1;}
      double PriceByBSFormula(double S0, double sigma, double r);

};

class EurPut
{
    public:
      double K, T; // strike
      int m;
      EurPut(double T_, double K_)
         {T=T_; K=K_; m=1;}
      double PriceByBSFormula(double S0, double sigma, double r);

};

class DifferenceOfOptions: public PathDepOption
{
   public:
      PathDepOption *ptr1,*ptr2;
      double Payoff(SamplePath &S)
      {
         return ptr1->Payoff(S)-ptr2->Payoff(S);
      }
      DifferenceOfOptions(double T_, int m_,
            PathDepOption *ptr1_,
            PathDepOption *ptr2_)
      {T=T_; m=m_; ptr1=ptr1_;ptr2=ptr2_;}
};
#endif
