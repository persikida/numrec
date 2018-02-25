#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <functional>
#include <array>
#include <algorithm>

const auto pi = 3.141592653589793;

using Arr3 = std::array<double, 3>;

std::vector<std::array<double, 3>> solve(
	const std::function<double(double)> &p,
	const std::function<double(double)> &q,
	const std::function<double(double)> &f,
	const double &t0,
	const double &tend,
	const Arr3 &bound_s,
	const Arr3 &bound_e
);

int main()
{

	auto p5 = [](const double &x)
	{
		return 0;
	};
	auto q5 = [](const double &x)
	{
		return -2.0 / (std::pow(std::cos(x), 2));
	};
	auto f5 = [](const double &x)
	{
		return 0;
	};
	std::ofstream out("5sh.csv");
	auto sol = solve(p5, q5, f5, 0, pi / 4, Arr3{ 0,1,1 }, Arr3{ 0,1,2 });
	for (std::size_t i = 0; i<sol.size(); i++)
		out << sol[i][0] << ' ' << sol[i][1] << ' ' << std::tan(sol[i][0]) << '\n';
	out.close();
	out.open("17_sh.csv");
	auto p17 = [](const double &x)
	{
		return (x - 3) / (x*x - 1);
	};
	auto q17 = [](const double &x)
	{
		return -1.0 / (x*x - 1);
	};
	auto f17 = [](const double &x)
	{
		return 0;
	};
	sol = solve(p17, q17, f17, 0, 1, Arr3{ 0,1,0 }, Arr3{ 1,1,-0.75 });
	for (std::size_t i = 0; i<sol.size(); i++)
	{
		auto x = sol[i][0];
		out << x << ' ' << sol[i][1] << ' ' << x - 3 + (1.0 / (x + 1)) << '\n';
	}
}

std::vector<std::array<double, 3>> solve(
	const std::function<double(double)> &p,
	const std::function<double(double)> &q,
	const std::function<double(double)> &f,
	const double &t0,
	const double &tend,
	const Arr3 &bound_s,
	const Arr3 &bound_e
)
{
	auto U_ = [](double x, double u, double v)
	{
		return v;
	};
	auto V_ = [&](double x, double u, double v)
	{
		return f(x) - p(x)*v - q(x)*u;
	};
	auto solveK = [&](
		const std::function<double(double, double, double)>&firstFunc,
		const std::function<double(double, double, double)>&secondFunc,
		const std::array<double, 3> &initCond,
		const double &tend,
		const double &h)
	{
		{
			std::vector<std::array<double, 3>> vec(int(std::abs((tend - initCond[0]) / h)) + 1);

			vec[0] = initCond;

			for (size_t count = 1; count<vec.size(); count++)
			{
				auto i = initCond[0] + h*count;
				auto k1 = h * firstFunc(i, vec[count - 1][1], vec[count - 1][2]);
				auto l1 = h * secondFunc(i, vec[count - 1][1], vec[count - 1][2]);
				auto k2 = h * firstFunc(i + h / 2, vec[count - 1][1] + k1 / 2, vec[count - 1][2] + l1 / 2);
				auto l2 = h * secondFunc(i + h / 2, vec[count - 1][1] + k1 / 2, vec[count - 1][2] + l1 / 2);
				auto k3 = h * firstFunc(i + h / 2, vec[count - 1][1] + k2 / 2, vec[count - 1][2] + l2 / 2);
				auto l3 = h * secondFunc(i + h / 2, vec[count - 1][1] + k2 / 2, vec[count - 1][2] + l2 / 2);
				auto k4 = h * firstFunc(i + h, vec[count - 1][1] + k3, vec[count - 1][2] + l3);
				auto l4 = h * secondFunc(i + h, vec[count - 1][1] + k3, vec[count - 1][2] + l3);

				vec[count][0] = i;
				vec[count][1] = vec[count - 1][1] + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
				vec[count][2] = vec[count - 1][2] + (l1 + 2 * l2 + 2 * l3 + l4) / 6;
			}
			return vec;
		}
	};

	double alpha1 = bound_s[0];
	double alpha2 = bound_e[0];
	double beta1 = bound_s[1];
	double beta2 = bound_e[1];
	double A = bound_s[2];
	double B = bound_e[2];

	auto tol = 2.220446049250313e-16;
	auto mode = std::abs(alpha1) <= tol && std::abs(A) <= tol;

	double nu = -1;
	double delta = 0.1;
	double prevNu = nu + delta;

	auto F = [&](double nu)
	{
		double v0, u0;
		if (mode)
		{
			v0 = (B - alpha2*nu) / beta2;
			u0 = nu;
			std::array<double, 3> temp = { tend,u0,v0 };
			auto value = solveK(U_, V_, temp, t0, -0.001).back();
			auto delta = alpha1*value[1] + beta1*value[2] - A;
			return delta;

		}
		else
		{
			v0 = (A - alpha1*nu) / beta1;
			u0 = nu;
			std::array<double, 3> temp = { t0,u0,v0 };
			auto value = solveK(U_, V_, temp, tend, 0.001).back();
			auto delta = alpha2*value[1] + beta2*value[2] - B;
			return delta;
		}

	};
	size_t count = 1;
	auto valF = F(nu);
	while (true)
	{

		auto temp = nu - valF*(nu - prevNu) / (valF - F(prevNu));
		prevNu = nu;
		nu = temp;

		valF = F(nu);
		if (abs(valF) <= 0.0001)
		{
			if (mode)
			{
				auto vec = solveK(U_, V_, { tend,nu,(B - alpha2*nu) / beta2 }, t0, -0.001);
				std::reverse(vec.begin(), vec.end());
				return vec;
			}
			else
				return solveK(U_, V_, { t0,nu,(A - alpha1*nu) / beta1 }, tend, 0.001);


		}
		if (++count > 1000)
		{
			std::cerr << "Calculation error\n";
			return std::vector<std::array<double, 3>>();
		}
	}

}



