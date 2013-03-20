/*
 * LBPTest.cpp
 * Example routines for using the LBP class.
 *
 *  Created on: Jan 25, 2013
 *      Author: Navid Nourani-Vatani
 *      Email: Navid.Nourani-Vatani@sydney.edu.au
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <opencv/cv.h>

#include "LBP.hpp"

using namespace lbp;

/**
 * Load and image and calculate the LBP-HF descriptor for the whole image
 */
void example_1( void ) {
	// Read an (RGB) image and convert to monochrome
	cv::Mat img = imread( "./test_image_1.bmp", 0 );
	// convert to double precision
	img.convertTo( img, CV_64F );

	// Create an LBP instance of type HF using 8 support points
	LBP lbp( 8, LBP_MAPPING_HF );
	// Calculate the descriptor
	lbp.calcLBP( img );
	// Calculate Fourier tranformed histogram
	vector<double> hist = lbp.calcHist().getHist( false );

	// Print out the histogram values
	double sum = 0;
	cout << "hist = [";
	for( int i = 0; i < hist.size(); i++ ) {
		cout << hist[i] << ", ";
		sum += hist[i];
	}
	cout << "]; " << endl;
	cout << "hist sum=" << sum << endl;
}

/**
 * Load an image, calculate LBP riu2 and then calculate the histogram on sub-regions of the image
 */
void example_2( void ) {
	// Read an (RGB) image and convert to monochrome
	cv::Mat img = imread( "/Users/navid/proj/LBP/test_image_1.bmp", 0 );
	// convert to double precision
	img.convertTo( img, CV_64F );
    

	// Create an LBP instance of type rotation invariant uniform 2 using 8 support points
	LBP lbp( 8, LBP_MAPPING_RIU2 );
	// Calculate the descriptor image and get it
	Mat lbpImg = lbp.calcLBP( img ).getLBPImage();

	// Get the histogram for sub-images
	for( int j = 0; j < 2; j++ ) {
		for( int i = 0; i < 2; i++ ) {
			// Get a sub-image the size of 1/4 of the whole image
			int x = lbpImg.cols / 2 * i;
			int y = lbpImg.rows / 2 * j;
			int w = lbpImg.cols/2-2;
			int h = lbpImg.rows/2-2;
			Mat subImg( lbpImg, Range(y,y+h), Range(x,x+w) );

			vector<double> hist = lbp.calcHist( &subImg ).getHist();
			// Print out the histogram values
			cout << "hist(" << j << "," << i << ") = [";
			for( int i = 0; i < hist.size(); i++ ) {
				cout << hist[i] << ", ";
			}
			cout << "]; " << endl;
		}
	}

}

/** 
 * Calculate a mapping, save it and load it. This is especially useful for 
 *	larger mappings (24 pts) which can takes many seconds to calculate.
 */
void example_3( void ) {
	clock_t startTime, endTime;
	
	LBP lbp( 16, LBP_MAPPING_U2 );
    cout << lbp.toString() << endl;
    startTime = clock();
    lbp.saveMapping( "mapping.txt" );
	endTime = clock();
	cout << "save took " << double( endTime - startTime ) / double( CLOCKS_PER_SEC ) << "s"
	<< endl;
	
    LBP lbp2;
    startTime = clock();
    lbp2.loadMapping("mapping.txt");
    endTime = clock();
    cout << lbp2.toString() << endl;
    cout << "load took " << double( endTime - startTime ) / double( CLOCKS_PER_SEC ) << "s"
    << endl;

	
}
int main( int argc, char ** argv ) {

	clock_t startTime, endTime;

	startTime = clock();
	example_1();
	endTime = clock();
	cout << "Example 2 took " << double( endTime - startTime ) / double( CLOCKS_PER_SEC ) << "s"
				<< endl;
    
    
	return 0;
}

