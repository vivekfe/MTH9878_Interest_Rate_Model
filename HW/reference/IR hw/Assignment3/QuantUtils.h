#include <cmath>
#include <random>

using namespace System;

/*
 * The option struct
 */
struct Option
{
	double T, F_0, K;
	bool isCall;
	// function to evaluate the option price
	double (*normalPricer)(double, double, double, double, bool);

	// function to evaluate the vega
	double (*pricer_prime)(double, double, double, double, bool);

	Option(double t, double f_0, double k, bool call, double (*normalP)(double, double, double, double, bool),
		double (*p_prime)(double, double, double, double, bool))
	{
		T = t;
		F_0 = f_0;
		K = k;

		isCall = call;

		normalPricer = normalP;
		pricer_prime = p_prime;
	}

	/*
	 * Calculate the option price using the given parameters
	 */
	double GetNormalOptionPriceWithVol(double vol)
	{
		return normalPricer(T, F_0, K, vol, isCall);
	}

	/*
	 * Calculate the vega using the given parameters
	 */
	double GetNormalOptionPriceWithVolDeriv(double vol)
	{
		return pricer_prime(T, F_0, K, vol, isCall);
	}
};

class QuanUtils
{
private:
	/*
	 * Use the Simpson's method to evaluate the cumulated normal distribution
	 */
	static double SimpsonApproximationCore(double a, double b, int n, double (*function)(double))
	{
		  double h = (b - a) / n;

		  double output = function(a) / 6 + function(b) / 6;

		  for (int i = 1; i < n; i++)
		  {
			output = output + function(a + i * h) / 3;
		  }

		  for (int i = 1; i <= n; i++)
		  {
			output = output + 2 * function(a + (i - 0.5) * h) / 3;
		  }

		  return h * output;
	}

public:
	/*
	 * Generate the standard normal random variable
	 */
	static double GetStandardNormal()
	{
		static std::random_device rd;
		static std::mt19937 eng(rd());
		static std::normal_distribution<> dist( 0.0, 1.0 ) ;

		return dist(eng) ;
	}

	/*
	 * Generate two correlated standard normal random variable
	 * 
	 * @param: rho: correlation between two normal
	 * @param: z1: the first standard normal variable
	 * @param: z2: the second standard normal variable
	 */
	static void GetCorrelatedStandardNormal(double rho, double& z1, double& z2)
	{
		z1 = GetStandardNormal();

		z2 = GetStandardNormal();

		z2 = rho * z1 + sqrt(1 - rho * rho) * z2;
	}

	/*
	 * Generate the density of the input standard normal variable
	 * 
	 * @param: z: the standard normal variable
	 */
	static double GetNormalDist(double z)
	{
        double result = 1.0 / Math::Sqrt(2.0 * Math::PI) * Math::Exp(-1.0 * Math::Pow(z, 2) / 2);

        return result;
	}

	/*
	 * return the cumulative normal value
	 * 
	 * @param: z: the standard normal variable
	 */
	static double CumulativeNormal(double z)
	{
		if (z == 0)
			return 0.5;

		int n = 4;

		double tol = Math::Pow(10, -12);

		double result = SimpsonApproximation(0, Math::Abs(z), n, GetNormalDist, tol);

		if (z > 0)
			return 0.5 + result;
		else
			return 0.5 - result;
	}

	/*
	 * return Simpson approximation of the integral
	 * 
	 * @param: a: The lower bound of the integral
	 * @param: b: The upper bound of the integral
	 * @param: n: The starting partition
	 * @param: function: the integrand of the integral
	 */
	static double SimpsonApproximation(double a, double b, int n, double (*function)(double), double tol)
	{
		 int finalFactor = 1;
		 if (b < a)
		 {
			   double c = a;

				a = b;

				b = c;

				finalFactor = -1;
		  }

		  double result = SimpsonApproximationCore(a, b, n, function);

		  double nextResult = result;

		  do
		  {
			   result = nextResult;

			   n = n * 2;

				nextResult = SimpsonApproximationCore(a, b, n, function);
		  }
		  while (Math::Abs(nextResult - result) >= tol);

		  return nextResult * finalFactor;
	}

	/*
	 * return impiled normal volatility using the option parameters
	 * 
	 * @param: x_0: Initial guess of the Newton method
	 * @param: option: The upper bound of the integral
	 * @param: targetVol: the target lognormal volality
	 * @param: tol: the tolarence
	 */
	static double GetNormalVolWithNewtonMethod(double x_0, Option &option, double targetVol, double tol)
    {
		double x_new = x_0;
        double x_old = x_0 - 1;

        while (Math::Abs(x_new - x_old) > tol)
        {
            x_old = x_new;

            x_new = x_old - (option.GetNormalOptionPriceWithVol(x_old) - targetVol) / option.GetNormalOptionPriceWithVolDeriv(x_old);
        }

        return x_new;
    }
};