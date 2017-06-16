#pragma once
#include <fitsio.h>
#include <iostream>
#include <boost/multi_array.hpp>

namespace Library
{
	template <typename T>
	class ImageFactory
	{
	public:
		static int readImage(std::string filename, boost::multi_array<T, 2>& data,
			bool transposeImageAfterReading = false);
		static void writeImage(std::string filename, boost::multi_array<T, 2>& data,
				std::size_t rows, std::size_t cols, bool overwrite = true,
				std::string importHeaderFrom = "", bool transposeImageBeforeWriting = false);
		static void getDimensions(std::string filename, std::size_t& rows, std::size_t& cols);
	};

	template <typename T>
	int ImageFactory<T>::readImage(std::string filename, boost::multi_array<T, 2>& data,
		bool transposeImageAfterReading)
	{
//			int rows = dims[0]; // = max(y);
//			int cols = dims[1]; // = max(x);
//			// x: [0, cols]; where cols is the number of columns in the image
//			// y: [0, rows]; where rows is the number of rows in the image
		fitsfile* fitsIn;
		std::size_t rows, cols, w = 0, dims[2] = {0, 0};
		float* buffer;
		int nullVal = 0, anyNull = 0, status = 0;
		if(fits_open_data(&fitsIn, filename.c_str(), READONLY, &status))
		{
			fits_report_error(stderr, status);
			return -1;
		}
		if(fits_get_img_size(fitsIn, 2, (long*)dims, &status))
		{
			fits_report_error(stderr, status);
			return -1;
		}

		rows = dims[1];
		cols = dims[0];
		if(transposeImageAfterReading)
			std::swap(rows, cols);

		buffer = new float[rows * cols];
		if(fits_read_img(fitsIn, TFLOAT, 1, rows * cols, &nullVal, buffer,
				&anyNull, &status))
		{
			fits_report_error(stderr, status);
			return -1;
		}
		if(fits_close_file(fitsIn, &status))
		{
			fits_report_error(stderr, status);
			return -1;
		}
		data.resize(boost::extents[rows][cols]);
		for(size_t col = 0; col < cols; col++)
			for(size_t row = 0; row < rows; row++, w++)
				data[row][col] = buffer[w];
		delete buffer;
		delete fitsIn;
		return rows;
	}

	template <typename T>
	void ImageFactory<T>::writeImage(std::string filename, boost::multi_array<T, 2>& data,
			std::size_t rows, std::size_t cols, bool overwrite,
			std::string importHeaderFrom, bool transposeImageBeforeWriting)
	{
		std::size_t naxes[2] = { cols, rows }, w = 0;
		int status = 0;
		float* buffer = new float[rows * cols];
		fitsfile* fitsOut;
//		fitsfile* fitsReference;

		if(transposeImageBeforeWriting)
		{
			std::swap(naxes[0], naxes[1]);
			for(size_t col = 0; col < cols; col++)
				for(size_t row = 0; row < rows; row++, w++)
					buffer[w] = data[row][col];
		}
		else
			for(size_t row = 0; row < rows; row++)
				for(size_t col = 0; col < cols; col++, w++)
					buffer[w] = data[row][col];

		if(overwrite)
			remove(filename.c_str());
		if(fits_create_file(&fitsOut, filename.c_str(), &status))
			fits_report_error(stderr, status);

/* *** NOT SUPPORTED YET *** */
//		if(importHeaderFrom != "")
//		{
//			fits_open_file(&fitsReference, importHeaderFrom.c_str(), READONLY, &status);
//
//			if(fits_copy_header(fitsReference, fitsOut, &status))
//				fits_report_error(stderr, status);
//			if(fits_close_file(fitsReference, &status))
//					fits_report_error(stderr, status);
//			delete fitsReference;
//		}

		if(fits_create_img(fitsOut, FLOAT_IMG, 2, (long*)naxes, &status))
			fits_report_error(stderr, status);
		if(fits_write_img(fitsOut, TFLOAT, 1, rows * cols, buffer, &status))
			fits_report_error(stderr, status);
		if(fits_close_file(fitsOut, &status))
			fits_report_error(stderr, status);
//		delete fitsOut;
	}

	template <typename T>
	void ImageFactory<T>::getDimensions(std::string filename, std::size_t& rows, std::size_t& cols)
	{
		fitsfile* fitsIn;
		int status = 0;
		size_t dims[2] = {0, 0};
		if(fits_open_data(&fitsIn, filename.c_str(), READONLY, &status))
			fits_report_error(stderr, status);
		if(fits_get_img_size(fitsIn, 2, (long*)dims, &status))
			fits_report_error(stderr, status);
		if(fits_close_file(fitsIn, &status))
			fits_report_error(stderr, status);
		delete fitsIn;
		cols = dims[0];
		rows = dims[1];
	}
}
