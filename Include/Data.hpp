/*
 * Data.hpp
 *
 *  Created on: 20 giu 2017
 *      Author: enzo
 */

#pragma once

#include "Maintainance.hpp"
#include <vector>

using namespace std;

namespace Library
{
	namespace Data
	{
		template <typename T>
		class Data
		{
		public:
			Data() { }
			Data(vector<T>& data, ushort_t nDataPoints);
			Data(ushort_t nDataPoints);
			Data(Data<T>* other);
			Data(Data<T>& other);
			T& operator[](ushort_t index) { return treatedData[index]; }
			vector<T>& getRawData() { return rawData; }
			vector<T>& getTreatedData() { return treatedData; }
			vector<T>& getDeviation() { return deviation; }
			ushort_t getNDataPoints() { return nDataPoints; }
			void initData();
		private:
			vector<T> rawData;
			vector<T> treatedData;
			vector<T> deviation;
			ushort_t nDataPoints;
			T zeroLevel;
		};

		template <typename T>
		Data<T>::Data(vector<T>& data, ushort_t nDataPoints) : rawData(nDataPoints),
			treatedData(nDataPoints), deviation(nDataPoints)
		{
			this->nDataPoints = nDataPoints;
			for(ushort_t x = 0; x < nDataPoints; x++)
				rawData[x] = data[x];
		}

		template <typename T>
		Data<T>::Data(ushort_t nDataPoints) : rawData(nDataPoints),
			treatedData(nDataPoints), deviation(nDataPoints)
		{
			this->nDataPoints = nDataPoints;
		}

		template <typename T>
		Data<T>::Data(Data<T>& other) : rawData(other.getNDataPoints()),
			treatedData(other.getNDataPoints()), deviation(other.getNDataPoints())
		{
			this->nDataPoints = other.getNDataPoints();
			for(ushort_t x = 0; x < nDataPoints; x++)
			{
				rawData[x] = other.getRawData()[x];
				treatedData[x] = other.getTreatedData()[x];
				deviation[x] = other.getDeviation()[x];
			}
		}

		template <typename T>
		Data<T>::Data(Data<T>* other) : rawData(other->getNDataPoints()),
			treatedData(other->getNDataPoints()), deviation(other->getNDataPoints())
		{
			this->nDataPoints = other->getNDataPoints();
			for(ushort_t x = 0; x < nDataPoints; x++)
			{
				rawData[x] = other->getRawData()[x];
				treatedData[x] = other->getTreatedData()[x];
				deviation[x] = other->getDeviation()[x];
			}
		}

		template <typename T>
		void Data<T>::initData()
		{
			T max = rawData[0], y = 0, y2 = 0, points = 0, gain = 0.68;
			int semi_width = 12;
			T center = 0;
			for(int i = 1; i < nDataPoints; i++)
			{
				if(rawData[i] > max)
				{
					max = rawData[i];
					center = i;
				}
			}
			for(int i = 0; i < nDataPoints; i++)
			{
				if(i < center - semi_width || i > center + semi_width)
				{
					y += rawData[i];
					y2 += rawData[i] * rawData[i];
					points++;
				}
			}
			T var = (y2 / points) - y * y / (points * points);
			for(int i = 0; i < nDataPoints; i++)
			{
				if(i > center - semi_width && i < center + semi_width)
				{
					this->treatedData[i] = rawData[i];
					this->deviation[i] = sqrt(fabs(rawData[i]) / gain + var);
				}
				else
				{
					this->treatedData[i] = sqrt(var);
					this->deviation[i] = sqrt(var);
				}
			}
			this->zeroLevel = sqrt(var);
		}
	}
}
