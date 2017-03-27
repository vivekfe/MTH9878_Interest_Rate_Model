//
//
//                                                                    Test.h
//

#ifndef TEST_H
#define TEST_H

#include "xlw/MyContainers.h"
#include <xlw/CellMatrix.h>
#include <xlw/DoubleOrNothing.h>
#include <xlw/ArgList.h>
#include <cmath>
#include "../QuantUtils.h"

using namespace xlw;



//<xlw:libraryname=MyInterestRateModelLibrary


double // Normal model option value
NormalOptionValue(double T // The maturity of the option
				, double F_0 // The current forward swap rate
				, double K   // The strike 
				, double sigma // The volatility of the underlying
				, bool isCall // Whether it is a call or put
				);

double // Calculate the vega of the normal model
NormalOptionValueDeriv(double T // The maturity of the option
					, double F_0 // The current forward swap rate
					, double K   // The strike 
					, double sigma // The volatility of the underlying
					, bool isCall // Whether it is a call or put
					);

double //Lognormal model option value
LognormalOptionValue(double T // The maturity of the option
					, double F_0 // The current forward swap rate
					, double K   // The strike 
					, double sigma // The volatility of the underlying
					, bool isCall // Whether it is a call or put
					);

double //Convert Lognormmal vol to Normal vol using the Newton's method
ConvertLognormalToNormalVol(double lognormalOptionPrice  // The option price from the Black model 
						, double T // The maturity of the option
						, double F_0 // The current forward swap rate
						, double K // The strike
						, bool isCall // Whether it is a call or put
						);

double //Get the SABR implied normal volatility
SABRImpliedNormalVol(double T // The maturity of the option
					, double F_0 // The current forward swap rate
					, double K // The strike
					, double sigma // The vol_0
					, double alpha // alpha in the SABR model
					, double beta // beta in the SABR model
					, double rho // rho in the SABR model
					);

double //Get the SABR implied log normal volatility
SABRImpliedLogNormalVol(double T // The maturity of the option
					, double F_0 // The current forward swap rate
					, double K // The strike
					, double sigma // The vol_0
					, double alpha // alpha in the SABR model
					, double beta // beta in the SABR model
					, double rho // rho in the SABR model
					);

double //Monte Carlo Simulation for SABR model
MonteCarloGetSABROptionPrice(double T
					, double F_0 // The current forward swap rate
					, double K // The strike
					, double sigma // The vol_0
					, double alpha // alpha in the SABR model
					, double beta // beta in the SABR model
					, double rho // rho in the SABR model
					, bool isCall //Whether it is a call or put
					, int numOfPath // Number of Monte Carlo simulation paths
					, int stepsInYear // Number of time step in one path
					);

#endif
