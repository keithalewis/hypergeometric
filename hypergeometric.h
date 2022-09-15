#pragma once
#include <cmath>
#include <limits>

namespace fms {

	// sum_0^infty (a)_n (b)_n/(c)_n x^n/n!
	// (q)_n = q(q + 1)...(q + n - 1), (q)_0 = 1
	template<class X = double>
	inline X _2F1(X a, X b, X c, X x)
	{
		X eps = sqrt(std::numeric_limits<X>::epsilon());
		X ab_cn = 1; // (a)_n (b)_n/(c)_n
		X x_n = 1;   // x^n/n!
		X s = 1;     // n = 0

		X ds = 1;
		for (int n = 1; fabs(ds) > eps || fabs(x_n) > eps; ++n, ++a, ++b, ++c) {
			X b_c = (b != 0 && c != 0) ? b / c : 1;
			ab_cn *= a * b_c;
			x_n *= x / n;
			ds = ab_cn * x_n;
			s += ds;
		}

		return s;
	}

#ifdef _DEBUG
	template<class X = double>
	inline bool test_2F1()
	{
		X eps = sqrt(std::numeric_limits<X>::epsilon());
		{
			// _2F1(1,b,b,x) = (1 - x)^-1 since (1)_n = n!
			X one = 1;
			X b = 1;
			X x = 0.1;
			X F = _2F1(one, b, b, x);
			X dF = F - 1 / (1 - x);
			if (fabs(dF) > eps) {
				return false;
			}
		}

		return true;
	}
#endif // _DEBUG
}
