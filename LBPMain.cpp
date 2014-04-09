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
				exit( 1 );
			}
		}
	}
	fileName = argv[argc - 1];

	// Read an (RGB) image and convert to monochrome
	Mat imgOrg = imread( fileName, -1 );
	// convert to double precision
	imgOrg.convertTo( imgOrg, CV_64F );

	Mat lbpImg;
	switch( imgOrg.channels() ) {
		case 1:
			lbpImg = Mat( imgOrg.size(), CV_8UC1, Scalar( 0 ) );
			break;
		case 3:
			lbpImg = Mat( imgOrg.size(), CV_8UC3, Scalar( 0 ) );
			break;
		default:
			cerr << "Unsupported number of image channels 1/3 only." << endl;
			exit( 1 );
	}

	// Create an LBP instance of type "mapping" using "pts" support points
	LBP lbp( pts, LBP::strToType( mapping ) );

	for( int i = 0; i < imgOrg.channels(); i++ ) {
		// Copy channel i
		Mat img( imgOrg.size(), imgOrg.depth(), 1 );
		const int from_to1[] = { i, 0 };
		mixChannels( &imgOrg, 1, &img, 1, from_to1, 1 );

		// Calculate the descriptor
		lbp.calcLBP( img, rad, true );

		// Copy lbp image
		const int from_to2[] = {0, i};
		Mat tmpImg = lbp.getLBPImage();
		mixChannels( &tmpImg, 1, &lbpImg, 1, from_to2, 1 );
	}

	if( strcmp( outFilename.c_str(), "" ) == 0 ) {
		char lbpType[1], lbpPts[2], lbpRad[2];
		sprintf( lbpType, "%d", lbp.getMapping() );
		sprintf( lbpPts, "%d", pts );
		sprintf( lbpRad, "%d", rad );

		outFilename = fileName.substr( 0, fileName.length() - 4 ) + "_LBP" + lbpType + "_" + lbpRad
		+ "_" + lbpPts + ".png";
	}

	if( outputHist ) {
		// Calculate Fourier tranformed histogram
		vector<double> hist = lbp.calcHist().getHist( normalizeHist );
		ofstream ofs;
		ofs.open( (outFilename.substr( 0, fileName.length() - 4 ) + ".txt").c_str(), ios::out );
		for( int i = 0; i < hist.size(); i++ ) {
			if( i > 0 ) ofs << ", ";
			ofs << hist[i];
		}
		ofs << endl;
		ofs.close();
	}
	else {
		imwrite( outFilename + ".png", lbpImg);
//		lbp.saveLBPImage( outFilename + ".png" );
	}

	return 0;
}
