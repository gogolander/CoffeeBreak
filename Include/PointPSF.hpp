/*
 * PointPSF.hpp
 *
 *  Created on: 20 mag 2017
 *      Author: Enzo
 */

#pragma once
#include "Maintainance.hpp"
#include "Point.hpp"
#include <string>

using namespace std;

namespace Library
{
	namespace Point
	{
		template <typename T>
		class PointPSF : public Point<T>
		{
		public:
			PointPSF()
			{
				this->point = vector<T>(6);
				this->nDimensions = 6;
			}
			PointPSF(PointPSF<T>& other);

			virtual ~PointPSF() { }
			virtual bool isValid();
			virtual void copyTo(Point<T>* target);
			virtual ushort_t stringToIndex(string term);
			virtual Point<T>* clone() { return new PointPSF<T>(this); }
			virtual Point<T>* newInstance() { return new PointPSF<T>(); }

			string getType() { return "Point.PSF"; }
			T fromFWHM(T width) { return width / GAUSS_FWHM; }
			T toFWHM(T dispersion) { return dispersion * GAUSS_FWHM; }
			T getFWHM() { return 2. * sqrt((pow(2, 1. / this->point["beta"]) - 1) / this->point["b"]); }
		};

		template <typename T>
		PointPSF<T>::PointPSF(PointPSF<T>& other)
		{
			for(ushort_t index = 0; index < this->nDimensions; index++)
				this->get(index) = other[index];
		}

		template <typename T>
		bool PointPSF<T>::isValid()
		{
			for(ushort_t i; i < this->nDimensions; i++)
				if(this->point[i] <= 0)
					return false;
			return true;
		}

		template <typename T>
		ushort_t PointPSF<T>::stringToIndex(string term)
		{
			transform(term.begin(), term.end(), term.begin(), ::tolower);
			ushort_t value = (term == "amplitude" || term == "amp" || term == "A" || term == "a") * 1 +
					(term == "b" || term == "B") * 2 +
					(term == "center" || term == "c" || term == "C") * 3 +
					(term == "beta" || term == "Beta") * 4 +
					(term == "sigma" || term == "dispersion") * 5 +
					(term == "zero" || term == "zeroLevel" || term == "zero level" ||
							term == "Zero" || term == "const") * 6;
			return value - 1;
		}

		template <typename T>
		void PointPSF<T>::copyTo(Point<T>* target)
		{
			if(target == NULL)
				target = new PointPSF<T>();
			for(ushort_t i = 0; i < this->nDimensions; i++)
				target->get(i) = this->point[i];
		}
	}
}
