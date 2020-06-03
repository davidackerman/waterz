#ifndef WATERZ_OPERATORS_H__
#define WATERZ_OPERATORS_H__

#include <functional>
#include <limits>
#include <cmath>
#include "MergeProviders.hpp"

template <typename ScoreFunction1, typename ScoreFunction2, template <typename> class Op>
class BinaryOperator : public ScoreFunction1, public ScoreFunction2 {

public:

	typedef typename MergeProviders<
			typename ScoreFunction1::StatisticsProviderType,
			typename ScoreFunction2::StatisticsProviderType>::Value
		StatisticsProviderType;

	typedef typename ScoreFunction1::ScoreType  ScoreType;

	template <typename RegionGraphType>
	BinaryOperator(
			RegionGraphType& regionGraph,
			const StatisticsProviderType& statisticsProvider) :
		ScoreFunction1(regionGraph, statisticsProvider),
		ScoreFunction2(regionGraph, statisticsProvider) {}

	template <typename EdgeIdType>
	inline ScoreType operator()(EdgeIdType e) {

		return _op(ScoreFunction1::operator()(e), ScoreFunction2::operator()(e));
	}

private:

	Op<ScoreType> _op;
};

template <typename ScoreFunction, template <typename> class Op>
class UnaryOperator : public ScoreFunction {

public:

	typedef typename ScoreFunction::StatisticsProviderType StatisticsProviderType;
	typedef typename ScoreFunction::ScoreType  ScoreType;

	template <typename RegionGraphType>
	UnaryOperator(
			RegionGraphType& regionGraph,
			const StatisticsProviderType& statisticsProvider) :
		ScoreFunction(regionGraph, statisticsProvider) {}

	template <typename EdgeIdType>
	inline ScoreType operator()(EdgeIdType e) {

		return _op(ScoreFunction::operator()(e));
	}

private:

	Op<ScoreType> _op;
};

template <typename T>
struct one_minus {
	T operator()(const T& x) const { return 1.0 - x; }
};
template <typename T>
using OneMinus = UnaryOperator<T, one_minus>;

template <typename T>
struct invert {
	T operator()(const T& x) const { return 1.0/x; }
};
template <typename T>
using Invert = UnaryOperator<T, invert>;

template <typename T>
struct square {
	T operator()(const T& x) const {

		return x*x;
	}
};
template <typename T>
using Square = UnaryOperator<T, square>;

template <typename T>
struct cosemSignedDistanceTransformOfRadius {
	T operator()(const T& x) const {
		float r = std::sqrt(x/3.14159265); //Radius (max distance to expect) assuming contact area is circular
		r = (r<35) ? r : 35;//DGA: For COSEM, 35 (sqrt(1225)) voxels is about the distance beyond which y = 128*tanh(d/50) + 127 saturates (here d is in nm). So we don't wan't to calculate anything beyond that
		r = 128 * std::tanh(r/12.5) + 127; 
		//std::cout<<"whattttttt"<<r<<std::endl;
		return r;
	}
};
template <typename T>
using CosemSignedDistanceTransformOfRadius = UnaryOperator<T, cosemSignedDistanceTransformOfRadius>;

template <typename T1, typename T2>
using Add = BinaryOperator<T1, T2, std::plus>;

template <typename T1, typename T2>
using Subtract = BinaryOperator<T1, T2, std::minus>;

template <typename T1, typename T2>
using Multiply = BinaryOperator<T1, T2, std::multiplies>;

template <typename T>
struct save_divide {

	T operator()(const T& a, const T& b) const {

		if (std::abs(b) <= std::numeric_limits<T>::min()) {

			if (std::signbit(a*b)) // a*b < 0
				return std::numeric_limits<T>::lowest();
			else
				return std::numeric_limits<T>::max();
		}

		return a/b;
	}
};

template <typename T1, typename T2>
using Divide = BinaryOperator<T1, T2, save_divide>;

template <typename T>
struct step {

	T operator()(const T& a, const T& b) const {

		if (a < b)
			return 0;
		return 1;
	}
};
template <typename T1, typename T2>
using Step = BinaryOperator<T1, T2, step>;

#endif // WATERZ_OPERATORS_H__

