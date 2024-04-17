#include <iostream>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <list>
#include <numeric>
#include <random>
#include <tuple>
#include <vector>

using vector = std::vector<double>;

const double X_MAX = 0;
const double X_MIN = std::_Pi;
const auto R = 5;
const auto M = (R - 1) / 2;
const auto K = 100;
const auto A = 0.5 / 2;
const vector H = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
const auto P = 0.95;
const auto EPSILON = 0.01;
const auto N =
std::ceil(std::log(1 - P) / std::log(1 - (EPSILON / (X_MIN - X_MAX))));

double random(const double l, const double r) 
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<> dis(l, r);
	return dis(gen);
}

double func(double x) 
{
	return sin(x) + 0.5;
}

double func_J(double h, double w, double d) 
{
	return h * w + (1 - h) * d;
}

vector generate_alpha(size_t r) 
{
	std::list<double> alpha;
	alpha.push_back(random(0, 1));

	for (size_t i = 1; i < (r - 1) / 2; ++i) 
	{
		double sum = std::accumulate(alpha.begin(), alpha.end(), 0.0);
		auto a = 0.5 * random(0, 1 - sum);
		alpha.push_back(a);
		alpha.push_front(a);
	}
	double sum = std::accumulate(alpha.begin(), alpha.end(), 0.0);
	auto a = 0.5 * (1 - sum);
	alpha.push_back(a);
	alpha.push_front(a);
	vector alpha_array(r);
	std::copy(alpha.begin(), alpha.end(), alpha_array.begin());
	return alpha_array;
}

double manhattan_distance_w(const vector& f_filter) 
{
	double s = 0;
	for (size_t k = 1; k < f_filter.size(); ++k)
		s += abs(f_filter[k] - f_filter[k - 1]);
	return s;
}

double manhattan_distance_d(const vector& f_filter, const vector& f_noise) 
{
	double s = 0;
	for (size_t k = 0; k < f_filter.size(); ++k)
		s += abs(f_filter[k] - f_noise[k]);
	return s;
}

double manhattan_distance(double w, double d) 
{
	return abs(w) + abs(d);
}

std::vector<double> noise_overlap(double a) 
{
	vector vec_noise;
	for (size_t k = 0; k < K; ++k) 
	{
		auto x_k = X_MAX + k * (X_MIN - X_MAX) / K;
		auto el = func(x_k) + random(-a, a);
		vec_noise.push_back(el);
	}
	return vec_noise;
}

vector filter(const vector& vec_noise, const vector& vec_alpha)
{
	vector func;
	for (size_t k = 0; k < K; ++k) 
	{
		auto j_start = k - M;
		auto j_end = k + M;
		if (j_start < 0 || j_end > K)
			continue;

		double p = 1;
		for (size_t j = j_start; j < j_end; ++j)
		{
			auto el = pow(vec_noise[j], vec_alpha[j + M + 1 - k]);
			p *= el;
		}
		func.push_back(p);
	}
	return func;
}

int main()
{
	srand(time(0));
	std::cout << std::fixed << std::setprecision(3) << std::right << std::right;
	auto f_n = noise_overlap(A);
	std::vector<std::pair<std::tuple<double, double, std::vector<double>>,
		std::tuple<double, double, std::vector<double>>>>
		vec_info;
	for (const auto& h : H) {
		vector min_alpha(3);
		vector min_f_f;
		auto min_w = 0.;
		auto min_d = 0.;
		auto min_J = 100.;
		auto min_dis_J = 100.;
		for (int i = 0; i < N; ++i)
		{
			auto alpha = generate_alpha(R);
			auto f_f = filter(f_n, alpha);
			auto w = manhattan_distance_w(f_f);
			auto d = manhattan_distance_d(f_f, f_n);
			auto J = func_J(h, w, d);
			auto dis_J = manhattan_distance(w, d);
			if (J < min_J)
			{
				min_alpha = alpha;
				min_dis_J = dis_J;
				min_J = J;
				min_w = w;
				min_d = d;
				min_f_f = f_f;
			}
		}
		vec_info.push_back(std::make_pair(std::make_tuple(min_J, min_dis_J, min_alpha), 
			std::make_tuple(min_w, min_d, min_f_f)));
	}

	for (size_t h = 0; h < vec_info.size(); ++h)
	{
		std::cout << "lambda =" << std::setw(6) << H[h] << std::setw(11)
			<< "min J =" << std::setw(7) << std::get<0>(vec_info[h].first)
			<< std::setw(15) << "min dis J =" << std::setw(7)
			<< std::get<1>(vec_info[h].first) << std::setw(15)
			<< "min alpha =" << std::setw(6)
			<< std::get<2>(vec_info[h].first)[0] << std::setw(6)
			<< std::setw(6) << std::get<2>(vec_info[h].first)[1]
			<< std::setw(6) << std::setw(6)
			<< std::get<2>(vec_info[h].first)[2] << std::setw(11)
			<< "min w =" << std::setw(7) << std::get<0>(vec_info[h].second)
			<< std::setw(11) << "min d =" << std::setw(7)
			<< std::get<1>(vec_info[h].second) << std::endl;
		std::cout << std::endl;
	}

	auto min = std::get<1>(vec_info[0].first);
	size_t min_h = 0;
	for (size_t h = 0; h < vec_info.size(); ++h)
	{
		min = std::get<1>(vec_info[h].first);
		if (std::get<1>(vec_info[h + 1].first) > min)
		{
			std::cout << "Optimum weight: " << H[h] << std::endl;
			std::cout << "Minimum distance: " << min << std::endl;
			min_h = H[h];
			break;
		}
	}
	system("pause");
	return 0;
}
