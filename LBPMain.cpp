/*
 * LBPMain.cpp
 *
 *  Created on: Aug 21, 2013
 *      Author: navid
 */

#include <cstdio>
#include <fstream>

#include "LBP.hpp"

using namespace std;
using namespace lbp;

void printHelp() {
	cout << "\nUsage: ./LBPMain [options] filename" << endl;
	cout << "\nOptions:" << endl;
	cout << "\t-r <int> - Radius (default=1)" << endl;
	cout << "\t-p <int> - Number of support points (default=8)" << endl;
	cout << "\t-m <string> - Mapping choose between: (default=none)" << endl;
	cout << "\t\tu2\n" << "\t\tri\n" << "\t\triu2\n" << "\t\thf" << endl;
	cout << "\t-o <string> - Output filename (default=filename_LBPm_r_p.png, "
				<< "where m,r,p correspond to mapping, radius and points.)" << endl;
	cout << "\t-h - Output histogram instead of LBP image" << endl;
	cout << "\t-hn - Output normalized histogram instead of LBP image" << endl;
}

int main( int argc, char ** argv ) {

	int rad = 1;
	int pts = 8;
	string mapping = "";
	string fileName = "";
	string outFilename = "";
	bool outputHist = false, normalizeHist = false;

	if( argc <= 2 ) {
		printHelp();
		exit( 1 );
	}
	else if( argc > 2 ) {
		// process arguments
		for( int i = 1; i < argc - 1; i++ ) {
			if( strcmp( argv[i], "-r" ) == 0 ) {
				rad = atoi( argv[i + 1] );
				i++;
			}
			else if( strcmp( argv[i], "-p" ) == 0 ) {
				pts = atoi( argv[i + 1] );
				i++;
			}
			else if( strcmp( argv[i], "-m" ) == 0 ) {
				mapping = argv[i + 1];
				i++;
			}
			else if( strcmp( argv[i], "-o" ) == 0 ) {
				outFilename = argv[i + 1];
				i++;
			}
			else if( strcmp( argv[i], "-h" ) == 0 ) {
				outputHist = true;
			}
			else if( strcmp( argv[i], "-hn" ) == 0 ) {
				outputHist = true;
				normalizeHist = true;
			}
			else {
				cerr << "invalid argument: \'" << argv[i] << "\'\n";
				printHelp();
				exit(1);
			}
		}
	}
	fileName = argv[argc - 1];

	// Read an (RGB) image and convert to monochrome
	cv::Mat img = imread( fileName, 0 );
	// convert to double precision
	img.convertTo( img, CV_64F );
	//cout << "image w/h = " << img.rows << "/" << img.cols << " (" << img.rows * img.cols << ")"
	//			<< endl;

	// Create an LBP instance of type HF using 8 support points
	LBP lbp( pts, lbp.strToType( mapping ) );
	// Calculate the descriptor
	lbp.calcLBP( img );

	if( strcmp( outFilename.c_str(), "" ) == 0 ) {
		char lbpType[1], lbpPts[2], lbpRad[2];
		sprintf( lbpType, "%d", lbp.getMapping() );
		sprintf( lbpPts, "%d", pts );
		sprintf( lbpRad, "%d", rad );

		outFilename = fileName.substr( 0, fileName.length() - 4 ) + "_LBP" + lbpType + "_" + lbpRad
					+ "_" + lbpPts;
	}

	if( outputHist ) {
		// Calculate Fourier tranformed histogram
		vector<double> hist = lbp.calcHist().getHist( normalizeHist );
		ofstream ofs;
		ofs.open( (outFilename + ".txt").c_str(), ios::out );
		for( int i = 0; i < hist.size(); i++ ) {
			if( i > 0 )	ofs << ", ";
			ofs << hist[i];
		}
		ofs << endl;
		ofs.close();
	}
	else {
		lbp.saveLBPImage( outFilename + ".png" );
	}

	return 0;
}
