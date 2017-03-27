
#include<cppinterface.h>

#pragma warning (disable : 4996)
using namespace System;

double // Normal model option value
NormalOptionValue(double T // The maturity of the option
				, double F_0 // The current forward swap rate
				, double K   // The strike 
				, double sigma // The volatility of the underlying
				, bool isCall // Whether it is a call or put
				)
{
	// Calculte the value of d+
	double d = (F_0 - K) / sigma / Math::Sqrt(T);

	// If it is a put, switch to d-
	if (!isCall)
		d *= -1.0;

	// Calculate the cumulative normmal distribution 
	double n = QuanUtils::CumulativeNormal(d);

	//Return the option price
	return sigma * Math::Sqrt(T) * (d * n + QuanUtils::GetNormalDist(d));
}

double // Calculate the vega of the normal model
NormalOptionValueDeriv(double T // The maturity of the option
					, double F_0 // The current forward swap rate
					, double K   // The strike 
					, double sigma // The volatility of the underlying
					, bool isCall // Whether it is a call or put
					)
{
	// Calculte the value of d+
	double d = (F_0 - K) / sigma / Math::Sqrt(T);

	if (!isCall)
		d *= -1.0;

	// Calculate the cumulative normmal distribution
	double n = QuanUtils::CumulativeNormal(d);

	double result = Math::Sqrt(T) * (d * n + QuanUtils::GetNormalDist(d));

	double d_prime = d / sigma * -1;

	result += d_prime * sigma * Math::Sqrt(T) * (n + d * QuanUtils::GetNormalDist(d) - d * QuanUtils::GetNormalDist(d));

	return result;
}

double //Lognormal model option value
LognormalOptionValue(double T // The maturity of the option
					, double F_0 // The current forward swap rate
					, double K   // The strike 
					, double sigma // The volatility of the underlying
					, bool isCall // Whether it is a call or put
					)
{
	//Calculate the value of d+
	double d_plus = (Math::Log(F_0 / K) + 0.5 * sigma * sigma * T) / sigma / Math::Sqrt(T);

	//Calculate the value of d-
	double d_minus = d_plus - sigma * Math::Sqrt(T);

	//Calculate the cumulative distribution of d+
	double n_plus = QuanUtils::CumulativeNormal(d_plus);

	//Calculate the cumulative distribution of d-
	double n_minus = QuanUtils::CumulativeNormal(d_minus);

	double result = 0;

	if (isCall)
	{
		result = F_0 * n_plus - K * n_minus;
	}
	else
	{
		result = K * (1 - n_minus) - F_0 * (1 - n_plus);
	}

	return  result;
}

double //Convert Lognormmal vol to Normal vol using the Newton's method
ConvertLognormalToNormalVol(double lognormalOptionPrice  // The option price from the Black model 
						, double T // The maturity of the option
						, double F_0 // The current forward swap rate
						, double K // The strike
						, bool isCall // Whether it is a call or put
						)
{
	//Initial guess
	double vol0 = 0.01;

	double tol = Math::Pow(10.0, -12);

	Option option(T, F_0, K, isCall, NormalOptionValue, NormalOptionValueDeriv);

	//Using the Newton's method to calculate the Implied Normal Vol
	double normalVol = QuanUtils::GetNormalVolWithNewtonMethod(vol0, option, lognormalOptionPrice, tol);

	return normalVol;
}

double //Get the SABR implied normal volatility
SABRImpliedNormalVol(double T // The maturity of the option
					, double F_0 // The current forward swap rate
					, double K // The strike
					, double sigma // The vol_0
					, double alpha // alpha in the SABR model
					, double beta // beta in the SABR model
					, double rho // rho in the SABR model
					)
{
	// Apply approximate formula in the percentage of moneyness is less than this value
	const double nearMoney = 0.05;

	double F_mid = (K + F_0) / 2;

	double gamma_1 = beta * Math::Pow(F_mid, beta - 1) / Math::Pow(F_mid, beta);

	double gamma_2 = beta * (beta - 1) * Math::Pow(F_mid, beta - 2) / Math::Pow(F_mid, beta);

	double C_F_mid = Math::Pow(F_mid, beta);

	double epsolon = T * alpha * alpha;

	// The extra term in the parenthesis
	double tmp = 1 + ((2 * gamma_2 - gamma_1 * gamma_1) /24
		* Math::Pow(sigma * C_F_mid / alpha, 2) + rho * gamma_1 / 4 * sigma * C_F_mid / alpha 
		+ (2 - 3 * rho * rho) / 24) * epsolon;

	double xi = alpha / sigma / (1 - beta) * (Math::Pow(F_0, 1.0 - beta) - Math::Pow(K, 1.0 - beta));

	double d_xi = Math::Log((Math::Sqrt(1 - 2 * rho * xi + xi * xi) + xi - rho) / (1 - rho));

	double result = 0;

	if (abs(F_0 - K) / K > nearMoney)
		result = alpha * (F_0 - K) / d_xi * tmp;
	else
		result = sigma * pow(F_0, beta) * tmp;

	return result;
}

double //Get the SABR implied log normal volatility
SABRImpliedLogNormalVol(double T // The maturity of the option
					, double F_0 // The current forward swap rate
					, double K // The strike
					, double sigma // The vol_0
					, double alpha // alpha in the SABR model
					, double beta // beta in the SABR model
					, double rho // rho in the SABR model
					)
{
	// Apply approximate formula in the percentage of moneyness is less than this value
	const double nearMoney = 0.05;

	double F_mid = (K + F_0) / 2;

	double gamma_1 = beta * Math::Pow(F_mid, beta - 1) / Math::Pow(F_mid, beta);

	double gamma_2 = beta * (beta - 1) * Math::Pow(F_mid, beta - 2) / Math::Pow(F_mid, beta);

	double C_F_mid = Math::Pow(F_mid, beta);

	double epsolon = T * alpha * alpha;

	// The extra term in the parenthesis
	double tmp = 1 + ((2 * gamma_2 - gamma_1 * gamma_1 + 1 / F_mid / F_mid) /24
		* Math::Pow(sigma * C_F_mid / alpha, 2) + rho * gamma_1 / 4 * sigma * C_F_mid / alpha 
		+ (2 - 3 * rho * rho) / 24) * epsolon;

	double xi = alpha / sigma / (1 - beta) * (Math::Pow(F_0, 1.0 - beta) - Math::Pow(K, 1.0 - beta));

	double d_xi = Math::Log((Math::Sqrt(1 - 2 * rho * xi + xi * xi) + xi - rho) / (1 - rho));

	double result = 0;

	if (abs(F_0 - K) / K > nearMoney)
		result = alpha * log(F_0 / K) / d_xi * tmp;
	else
		result = sigma * pow(F_0, beta - 1) * tmp;

	return result;
}

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
					)
{
	//const int stepsInYear = 360;

	double timeStep = 1.0 / stepsInYear;

	//Maturity of the option in days
	int mat = stepsInYear * T;

	double price = 0;

	for (int i = 0; i < numOfPath; i++)
	{
		// Beginning value of F and sigma
		double F_i = F_0;
		double sigma_i = sigma;

		int currTime = 1;

		while(currTime <= mat)
		{
			double x, y;

			//Generate two standard normal random variable with given correlation
			QuanUtils::GetCorrelatedStandardNormal(rho, x, y);

			// The change of F
			double tmp = sigma_i * pow(F_i, beta) * sqrt(timeStep) * x;

			F_i = Math::Max(F_i + tmp, 0.0);

			sigma_i = sigma_i * exp(alpha * sqrt(timeStep) * y - alpha * alpha * timeStep / 2);

			++currTime;
		}

		// Calculate the option payoff
		if (isCall)
		{
			price += Math::Max(F_i - K, 0.0);
		}
		else
			price += Math::Max(K - F_i, 0.0);
	}

	return price / numOfPath;
}

