/**
 * LBP.hpp
 * Implements the local binary pattern (LBP) texture descriptors
 *
 *  Created on: Jan 25, 2013
 *      Author: Navid Nourani-Vatani
 *      Email: Navid.Nourani-Vatani@sydney.edu.au
 *
 *  The methods implemented here are inspired by the Matlab code available
 *  from web site of the University of Oulu:
 *  	http://www.cse.oulu.fi/CMV/Downloads/LBPMatlab
 *  You should cite the appropriate publications when using this code.
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
 *
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <climits>
#include <cmath>
#include <complex>
#include <string>

#include <fftw3.h>

#include <opencv2/opencv.hpp>
#include <opencv2/highgui/highgui.hpp>


#ifndef _LBP_H_
#define _LBP_H_

#ifndef NULL
#define NULL 0
#endif

// SWIG doesn't understand 'using' directives properly
// so disable them if doing the SWIG pass
#ifndef SWIG
using namespace std;
using namespace cv;
#endif

// enable/disable use of mixed OpenCV API in the code below.
#define DEMO_MIXED_API_USE 0
#define SSTR( x ) dynamic_cast< std::ostringstream & >( \
( std::ostringstream() << std::dec << x ) ).str()

namespace lbp {
    
	enum MappingType {
		LBP_MAPPING_NONE = 0,
		LBP_MAPPING_U2,
		LBP_MAPPING_RI,
		LBP_MAPPING_RIU2,
		LBP_MAPPING_HF
	};
    
	static const string MappingTypeStr[] = {
				"none", "u2", "ri", "riu2", "hf" };
    
	class LBP {
        
	public:
		LBP();
		LBP( unsigned int samples, MappingType type );
		~LBP( void );
        
		/**
		 * Mapping methods
		 */
		LBP & generateMapping();
		LBP & generateMapping( unsigned int samples, MappingType type );

		bool saveMapping( string fileName );
		bool loadMapping( string fileName );
        
        static MappingType strToType( string s ) {
            
            if( s.compare( "u2" ) == 0 )
                return LBP_MAPPING_U2;
            else if( s.compare("ri") == 0 )
                return LBP_MAPPING_RI;
            else if ( s.compare("riu2") == 0 )
                return LBP_MAPPING_RIU2;
            else if( s.compare("hf")  == 0 )
                return LBP_MAPPING_HF;
            else
                return LBP_MAPPING_NONE;
        }
        
        MappingType getMapping(void) const {
        	return type;
        }

		/**
		 * Descriptor methods
		 */
		LBP & calcLBP( Mat img, double radius = 1., bool borderCopy=false );
		Mat getLBPImage( void ) const {
			return lbpImage;
		}
        
		bool saveLBPImage( string fileName );
		/**
		 * Histogram methods
		 */
		LBP & calcHist(void);
		LBP & calcHist( Mat mask );

		vector<double> getHist( bool norm = true );
		vector<double> constructHF( vector<double> h );
        
        
		/**
		 * Other methods
		 */
		std::string toString( bool verbose=false ) const;
        
	private:
        // Mapping variables
        MappingType type;
        vector<int> table;
        unsigned int samples;
        unsigned int num;
        // Fourier Histogram variables
        vector< vector<int> > orbits;
        double *fftIn;
        complex<double> *fftOut;
        fftw_plan fftPlan;
        unsigned int fftN;
        unsigned int fftHermN;
        vector<double> hf;
        // Histogram
        vector<double> h;
        // Descriptor variables
        Mat lbpImage;
        MatND hist;
        
        // Private bit operation methods
        int NumberOfSetBits( int i ) {
            i = i - ((i >> 1) & 0x55555555);
            i = (i & 0x33333333) + ((i >> 2) & 0x33333333);
            return (((i + (i >> 4)) & 0x0F0F0F0F) * 0x01010101) >> 24;
        }
        unsigned int rotateLeft( unsigned int i, unsigned int samples ) {
            unsigned int bg = ((i & (1 << (samples - 1))) >> (samples - 1)); // bitget(r,samples)
            unsigned int bs = (i << 1) & ((int) pow( 2., (int) samples ) - 1); // bitshift(r, 1, samples)
            unsigned int j = (bs + bg) & ((int) pow( 2., (int) samples ) - 1); // bitset( bs, 1, bg )
            return j;
        }
        int trailingZeroInd( unsigned int v ) {  // find the number of trailing zeros in v
            static const int Mod37BitPosition[] =  // map a bit value mod 37 to its position
            { 32, 0, 1, 26, 2, 23, 27, 0, 3, 16, 24, 30, 28, 11, 0, 13, 4, 7, 17, 0, 25,
                22, 31, 15, 29, 10, 12, 6, 0, 21, 14, 9, 5, 20, 8, 19, 18 };
            return Mod37BitPosition[(-v & v) % 37];
        }
        
        void initHF(void);

        LBP & calcHist( Mat * img, Mat * mask=NULL );
	};
    
}
;

#endif
