#include "stdafx.h"
#include "stb_image.h"
#include "stb_image_write.h"
#include <time.h>
#include <iostream>
#include <time.h>
#include <math.h>
#include <string.h>

#define QX_DEF_CHAR_MAX 255


#include <iomanip>

using namespace std;

/**
 * \brief This files gives an implementation of the unit test worker (C++ class Worker) for the HARDSI version of the Recursive Bilateral Filter.
 * \detail This code is based on an original implementation described in the paper from Qingxiong Yang: "Recursive Bilateral Filtering", European Conference on Computer Vision (ECCV) 2012, 399-413.
 */



// -------------------------------------------
// Global parameters
// -------------------------------------------
MATRIX_DEFINE(unsigned char, img_in);
MATRIX_DEFINE(unsigned char, IMG);

const float sigma_spatial = 0.12f;
const float sigma_range = 0.09f;


// --------------------------------------
// Auxiliary functions (header)
// --------------------------------------
inline void _recursive_bf(	float sigma_spatial, float sigma_range, int width, int height, int channel, float * buffer = 0)
{
	const int width_height			= width * height;
	const int width_channel			= width * channel;
	const int width_height_channel	= width * height * channel;

	rgb_t pixel, pixel_py, pixel_cy, pixel_x, texture_x;

	bool is_buffer_internal = (buffer == 0);
	if (is_buffer_internal)
		buffer = new float[(width_height_channel + width_height + width_channel + width) * 2];

	float * img_out_f		= buffer;
	float * img_temp		= &img_out_f[width_height_channel];
	float * map_factor_a	= &img_temp[width_height_channel];
	float * map_factor_b	= &map_factor_a[width_height]; 
	float * slice_factor_a	= &map_factor_b[width_height];
	float * slice_factor_b	= &slice_factor_a[width_channel];
	float * line_factor_a	= &slice_factor_b[width_channel];
	float * line_factor_b	= &line_factor_a[width];

	//compute a lookup table
	float range_table[QX_DEF_CHAR_MAX + 1];
	float inv_sigma_range = 1.0f / (sigma_range * QX_DEF_CHAR_MAX);
	for (int i = 0; i <= QX_DEF_CHAR_MAX; i++)
		range_table[i] = static_cast<float>(exp(-i * inv_sigma_range));

	float alpha = static_cast<float>(exp(-sqrt(2.0) / (sigma_spatial * width)));
	float ypr, ypg, ypb, ycr, ycg, ycb;
	float fp, fc;
	float inv_alpha_ = 1 - alpha;
	for (int y = 0; y < height; y++)
	{
		float			* temp_x	= &img_temp[y * width_channel];

		*temp_x++ = ypr = MATRIX_GET(IMG, 0,	y); 
		*temp_x++ = ypg = MATRIX_GET(IMG, 1,	y);
		*temp_x++ = ypb = MATRIX_GET(IMG, 2,	y);

		unsigned char tpr = MATRIX_GET(IMG, 0,	y);
		unsigned char tpg = MATRIX_GET(IMG, 1,	y);
		unsigned char tpb = MATRIX_GET(IMG, 2,	y);

		float * temp_factor_x = &map_factor_a[y * width];
		*temp_factor_x++ = fp = 1;

		// from left to right
		for (int x = 1; x < width; x++) 
		{
			unsigned char tcr = MATRIX_GET(IMG, x*channel,		y);
			unsigned char tcg = MATRIX_GET(IMG, x*channel+1,	y);
			unsigned char tcb = MATRIX_GET(IMG, x*channel+2,	y);

			unsigned char dr = abs(tcr - tpr);
			unsigned char dg = abs(tcg - tpg);
			unsigned char db = abs(tcb - tpb);
			int range_dist = (((dr << 1) + dg + db) >> 2);
			float weight = range_table[range_dist];
			float alpha_ = weight*alpha;

			unsigned char r = MATRIX_GET(IMG, x*channel,	y);
			unsigned char g = MATRIX_GET(IMG, x*channel+1,	y);
			unsigned char b = MATRIX_GET(IMG, x*channel+2,	y);

			*temp_x++ = ycr = inv_alpha_*r + alpha_*ypr; 
			*temp_x++ = ycg = inv_alpha_*g + alpha_*ypg; 
			*temp_x++ = ycb = inv_alpha_*b + alpha_*ypb;
			tpr = tcr; tpg = tcg; tpb = tcb;
			ypr = ycr; ypg = ycg; ypb = ycb;
			*temp_factor_x++ = fc = inv_alpha_ + alpha_*fp;
			fp = fc;
		}

		unsigned char r = MATRIX_GET(IMG, (width-1)*channel,	y);
		unsigned char g = MATRIX_GET(IMG, (width-1)*channel+1,	y);
		unsigned char b = MATRIX_GET(IMG, (width-1)*channel+2,	y);

		*--temp_x; *temp_x = 0.5f*((*temp_x) + b);
		*--temp_x; *temp_x = 0.5f*((*temp_x) + g);
		*--temp_x; *temp_x = 0.5f*((*temp_x) + r);

		tpr = b; 			tpg = g;			tpb = r;
		ypr = r;			ypg = r;			ypb = r;

		*--temp_factor_x; *temp_factor_x = 0.5f*((*temp_factor_x) + 1);
		fp = 1;

		// from right to left
		for (int x = width - 2; x >= 0; x--) 
		{
			unsigned char tcr = MATRIX_GET(IMG, x*channel,		y);
			unsigned char tcg = MATRIX_GET(IMG, x*channel+1,	y);
			unsigned char tcb = MATRIX_GET(IMG, x*channel+2,	y);

			unsigned char dr = abs(tcr - tpr);
			unsigned char dg = abs(tcg - tpg);
			unsigned char db = abs(tcb - tpb);
			int range_dist = (((dr << 1) + dg + db) >> 2);
			float weight = range_table[range_dist];
			float alpha_ = weight * alpha;

			unsigned char r = MATRIX_GET(IMG, x*channel,	y);
			unsigned char g = MATRIX_GET(IMG, x*channel+1,	y);
			unsigned char b = MATRIX_GET(IMG, x*channel+2,	y);

			ycr = inv_alpha_ * r + alpha_ * ypr; 
			ycg = inv_alpha_ * g + alpha_ * ypg; 
			ycb = inv_alpha_ * b + alpha_ * ypb;

			*--temp_x; *temp_x = 0.5f*((*temp_x) + ycr);
			*--temp_x; *temp_x = 0.5f*((*temp_x) + ycg);
			*--temp_x; *temp_x = 0.5f*((*temp_x) + ycb);
			tpr = tcr; tpg = tcg; tpb = tcb;
			ypr = ycr; ypg = ycg; ypb = ycb;

			fc = inv_alpha_ + alpha_*fp;
			*--temp_factor_x; 
			*temp_factor_x = 0.5f*((*temp_factor_x) + fc);
			fp = fc;
		}
	}
	alpha = static_cast<float>(exp(-sqrt(2.0) / (sigma_spatial * height)));
	inv_alpha_ = 1 - alpha;
	float * ycy, * ypy, * xcy;
	unsigned char * tcy, * tpy;
	memcpy(img_out_f, img_temp, sizeof(float)* width_channel);

	float * in_factor = map_factor_a;
	float*ycf, *ypf, *xcf;
	memcpy(map_factor_b, in_factor, sizeof(float) * width);
	for (int y = 1; y < height; y++)
	{
		xcy = &img_temp[y * width_channel];
		ypy = &img_out_f[(y - 1) * width_channel];
		ycy = &img_out_f[y * width_channel];

		xcf = &in_factor[y * width];
		ypf = &map_factor_b[(y - 1) * width];
		ycf = &map_factor_b[y * width];
		for (int x = 0; x < width; x++)
		{
			unsigned char pr = MATRIX_GET(IMG, x*channel,	y-1);
			unsigned char pg = MATRIX_GET(IMG, x*channel+1,	y-1);
			unsigned char pb = MATRIX_GET(IMG, x*channel+2,	y-1);

			unsigned char tr = MATRIX_GET(IMG, x*channel,	y);
			unsigned char tg = MATRIX_GET(IMG, x*channel+1,	y);
			unsigned char tb = MATRIX_GET(IMG, x*channel+2,	y);

			unsigned char dr = abs(tr - pr);
			unsigned char dg = abs(tg - pg);
			unsigned char db = abs(tb - pb);
			int range_dist = (((dr << 1) + dg + db) >> 2);
			float weight = range_table[range_dist];
			float alpha_ = weight*alpha;
			for (int c = 0; c < channel; c++) 
				*ycy++ = inv_alpha_*(*xcy++) + alpha_*(*ypy++);
			*ycf++ = inv_alpha_*(*xcf++) + alpha_*(*ypf++);
		}
	}
	int h1 = height - 1;
	ycf = line_factor_a;
	ypf = line_factor_b;
	memcpy(ypf, &in_factor[h1 * width], sizeof(float) * width);
	for (int x = 0; x < width; x++) 
		map_factor_b[h1 * width + x] = 0.5f*(map_factor_b[h1 * width + x] + ypf[x]);

	ycy = slice_factor_a;
	ypy = slice_factor_b;
	memcpy(ypy, &img_temp[h1 * width_channel], sizeof(float)* width_channel);
	int k = 0; 
	for (int x = 0; x < width; x++) {
		for (int c = 0; c < channel; c++) {
			int idx = (h1 * width + x) * channel + c;
			img_out_f[idx] = 0.5f*(img_out_f[idx] + ypy[k++]) / map_factor_b[h1 * width + x];
		}
	}

	for (int y = h1 - 1; y >= 0; y--)
	{
		xcy = &img_temp[y * width_channel];
		float*ycy_ = ycy;
		float*ypy_ = ypy;
		float*out_ = &img_out_f[y * width_channel];

		xcf = &in_factor[y * width];
		float*ycf_ = ycf;
		float*ypf_ = ypf;
		float*factor_ = &map_factor_b[y * width];
		for (int x = 0; x < width; x++)
		{
			unsigned char pr = MATRIX_GET(IMG, x*channel,	y+1);
			unsigned char pg = MATRIX_GET(IMG, x*channel+1,	y+1);
			unsigned char pb = MATRIX_GET(IMG, x*channel+2,	y+1);

			unsigned char tr = MATRIX_GET(IMG, x*channel,	y);
			unsigned char tg = MATRIX_GET(IMG, x*channel+1,	y);
			unsigned char tb = MATRIX_GET(IMG, x*channel+2,	y);

			unsigned char dr = abs(tr - pr);
			unsigned char dg = abs(tg - pg);
			unsigned char db = abs(tb - pb);
			int range_dist = (((dr << 1) + dg + db) >> 2);
			float weight = range_table[range_dist];
			float alpha_ = weight*alpha;

			float fcc = inv_alpha_*(*xcf++) + alpha_*(*ypf_++);
			*ycf_++ = fcc;
			*factor_ = 0.5f * (*factor_ + fcc);

			for (int c = 0; c < channel; c++)
			{
				float ycc = inv_alpha_*(*xcy++) + alpha_*(*ypy_++);
				*ycy_++ = ycc;
				*out_ = 0.5f * (*out_ + ycc) / (*factor_);
				*out_++;
			}
			*factor_++;
		}
		memcpy(ypy, ycy, sizeof(float) * width_channel);
		memcpy(ypf, ycf, sizeof(float) * width);
	}

	for (int y = 0; y < height; ++y)
	{
		for (int x = 0; x < width*channel; ++x)
		{
			MATRIX_SET(IMG, x, y, static_cast<unsigned char>(*img_out_f));
			img_out_f ++;
		}
	}

	if (is_buffer_internal)
		delete[] buffer;
}


inline void recursive_bf(	float sigma_spatial, float sigma_range,
							int width, int height, int channel,
							float * buffer = 0)
{
	_recursive_bf(sigma_spatial, sigma_range, width, height, channel, buffer);
}





// path where files are located, you may need to change this
const char images_folder_path[] = "./images/";

// test images:
const char file_name_testGirl[]		= "testGirl.jpg";		// size: 448 x 626
const char file_name_house[]		= "Thefarmhouse.jpg";	// size: 1440 x 1080
const char file_name_testpattern[]	= "testpatern5.png";	// size: 1920 x 1080


// timer uses 'test_runs' as divisor
class TestRunTimer 
{
	clock_t begTime;

public:
	void start() { begTime = clock(); }
	float elapsedTimeMS() { return float(clock() - begTime); }
};

// utility for setting output file name
template <size_t _Size>
char* modifyFilePath(char (&file_path)[_Size], const char* suffix)
{
	size_t l = strlen(file_path);
	// get rid of old extension
	for (size_t i = l - 1; i > 0; i--)
	{
		if (file_path[i] == '.')
		{
			file_path[i] = 0;
			break;
		}
	}

	// add current sigma values just for clarity
	char extra_text[64];
	sprintf(extra_text, "%0.3f_%0.3f", sigma_spatial, sigma_range);

	// add suffix
	strcat(file_path, "_");
	strcat(file_path, suffix);
	strcat(file_path, "_");
	strcat(file_path, extra_text);
	strcat(file_path, ".png"); // force PNG format

	return file_path;
}

// using original implementation, source code from
// https://github.com/ufoym/RecursiveBF
void testRunRecursiveBF_Original(const char* image_name)
{
	cout << "\nImage: " << image_name;
	char file_path[256];
	strcpy(file_path, images_folder_path);
	strcat(file_path, image_name);

	int width, height, channel;
	// Allocate and copy img_in
	unsigned char *img_in_no = stbi_load(file_path, &width, &height, &channel, 3);
	if (!img_in_no)
	{
		cout << "\nFailed to load image path: " << file_path;
		return;
	}
	cout << ", size: " << width << " x " << height;
	channel = 3; // require 3 channel for this test
	TestRunTimer timer;

	// memory reserve for filter algorithm before timer start
	float * buffer = new float[(width * height* channel + width * height + width * channel + width) * 2];

	// Allocate IMG and copy img_in
	MATRIX_ALLOCATE(unsigned char, width*channel, height, img_in);
	MATRIX_ALLOCATE(unsigned char, width*channel, height, IMG);
	for (int y = 0; y < height; ++y)
	{
		for (int x = 0; x < width*channel; ++x)
		{
			MATRIX_SET(img_in,	x, y, *img_in_no);
			MATRIX_SET(IMG,		x, y, *img_in_no++);
		}
	}

	timer.start();

	recursive_bf(	sigma_spatial, sigma_range,
					width, height,
					channel, buffer);

	cout << ", time ms: " << timer.elapsedTimeMS();
	
	delete[] buffer;

	modifyFilePath(file_path, "RBF");

	// Transform IMG to give it to stbi_write_png
	stbi_write_png(file_path, width, height, channel, IMG, width * 3);

	delete[] img_in;
	delete[] IMG;

	img_in	= NULL;
	IMG		= NULL;
}

/////////////////////////////////////////////////////////////////////////////

int main()
{
	cout << "test run \n";
	cout << fixed << setprecision(1);

	////////////////////////
	cout << "\nOriginal Recursive Bilateral Filter implementation";
	// image: testpattern
	testRunRecursiveBF_Original(file_name_testpattern);
	// image: house
	testRunRecursiveBF_Original(file_name_house);
	// image: testGirl
	testRunRecursiveBF_Original(file_name_testGirl);

	

	cout << "\nFinish";
	cin.get();

	return 0;
}

