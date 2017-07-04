/*
 * Mantainance.hpp
 *
 *  Created on: 30 mag 2017
 *      Author: Enzo
 */

#pragma once

#include <boost/multi_array.hpp>
namespace Library
{
	typedef unsigned short ushort_t;
	typedef unsigned short Index;
	template <typename T>
	using Vector = std::vector<T>;
	const double GAUSS_FWHM = 2.354820046; // 2 * sqrt(2 * log(2))
}
