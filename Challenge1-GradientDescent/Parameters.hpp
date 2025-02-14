#include "Grad.hpp"
#include <functional>
#include <vector>
#ifndef Parameters_library
#define Parameters_library

enum class method {Armijo,heavyBall,Nesterov}; //definition of all possible methods ot compute the minimum

class Parameters{ // definition of the class conteining all the parameters of the problem
    private:
        Point x0;  //inital guess of the minimum
        double tollr;  //tollerance over the residual
        double tolls;  //tollerance over the value of the gradient
        int maxIter;   //maximum number of iterations
        double alpha0; //value of alpha0
        double omega;  //value of omega for the condition in Armijo rule
        double mu;     //value of mu in the inverse decay rule
        double nu;     //value of nu in Nesterov aand heavy-ball methods
        std::function<double(Point)> f; //the function that needs to be minimized
        Grad G;   //the gradient of that function

    public:
        Parameters(Point X0,double Tollr,double Tolls,int MaxIter,double Alpha0,double Omega,double Mu,double Nu,std::function<double(Point)> F,Grad g):
        x0(X0),tollr(Tollr),tolls(Tolls),maxIter(MaxIter),alpha0(Alpha0),omega(Omega),mu(Mu),nu(Nu),f(F),G(g){}; //the constructor
        
        void print() const;  //displaies all the parameters

        Point getX0()const{return x0;};   //all getter needed inthe implementation to access to the private members
        double getTollr()const{return tollr;};
        double getTolls()const{return tolls;};
        int getMaxIter()const{return maxIter;};
        double getAlpha0()const{return alpha0;};
        double getOmega()const{return omega;};
        double getMu()const{return mu;};
        double getNu()const{return nu;};
        std::function<double(Point)> getF()const{return f;};
        Grad getG()const{return G;};

};


bool condition(const Point& x, double alpha,double omega,std::function<double(Point)> f,Grad G);  //checks if the condition of Armijo is satisfied

template <method M>
Point Minimum(const Parameters& p);  //main function that computes the minimum

#endif