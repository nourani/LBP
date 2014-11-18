#include <stdlib.h>

//#include <npp.h>
//#include <nppi.h>
#include <cuda.h>
#include <device_launch_parameters.h>
#include <cuda_runtime.h>

#define cimg_display 1
#include <CImg.h>

#include "LBPGPU.cuh"

__global__ void calcLBPKernel( const unsigned char * pSrc, unsigned char * pDst, const int width,
			const int height, const LBPMapping * mapping ) {
	// Get the index of the current core
	int tidX = blockIdx.x * blockDim.x + threadIdx.x;
	int tidY = blockIdx.y * blockDim.y + threadIdx.y;
	int idx = tidY * width + tidX;
	int idxX = idx / width;
	int idxY = idx % width;

	// Only continue if a distance 'radius' from the edge
	float r = mapping->radius;
	if( idxX < r || idxX >= width - r || idxY < r || idxY >= height - r )
		return;
//	if( idx != 7 )
//		return;

	const int samples = 8; //mapping->samples;
	float spoints[samples][2];
	float a = 2.f * M_PI / samples;
	for( int i = 0; i < samples; i++ ) {
		spoints[i][0] = +r * cos( float( i * a ) );
		spoints[i][1] = -r * sin( float( i * a ) );
	}

	float val=0;
	int idxC = idx;
	float valC = pSrc[idxC];
	for( int i = 0; i < samples; i++ ) {
		float x = spoints[i][0];
		float y = spoints[i][1];
		int ry = round( y );
		int rx = round( x );
		int indV = idxC + width * ry + rx;
		float valV = pSrc[indV];

		// Check if interpolation is needed.
		if( (fabs( x - rx ) > 1e-6) || (fabs( y - ry ) > 1e-6) ) {
			int fy = floor( y );
			int cy = ceil( y );
			int fx = floor( x );
			int cx = ceil( x );

			int idxV1 = idxC + width * fy + fx;
			int idxV2 = idxC + width * fy + cx;
			int idxV3 = idxC + width * cy + fx;
			int idxV4 = idxC + width * cy + cx;
			// Calculate the interpolation weights.
			float tx = x - fx;
			float ty = y - fy;
			float w1 = (1 - tx) * (1 - ty);
			float w2 = tx * (1 - ty);
			float w3 = (1 - tx) * ty;
			float w4 = tx * ty;
			valV = pSrc[idxV1]*w1 + pSrc[idxV2]*w2 + pSrc[idxV3]*w3 + pSrc[idxV4]*w4;
		}
		val += valV - valC >= 0 ? pow( 2.f, (float)i ) : 0;
	}

	pDst[idx] = (unsigned char) val;
//	printf("%d (%d,%d): %f\n", idx, idxX, idxY, val);
}

/**
 *  The Wrapper function
 */
void calcLBPGPU( const unsigned char * h_src, unsigned char * h_dst, const int width, const int height,
			const LBPMapping * mapping ) {

	unsigned char *d_Src = NULL, *d_Dst = NULL;

	cudaMalloc( &d_Src, sizeof(char) * height * width );
	cudaMalloc( &d_Dst, sizeof(char) * height * width );
	cudaMemset( (void *) d_Dst, 0, sizeof(char) * height * width );
	cudaMemcpy( d_Src, h_src, sizeof(char) * height * width, cudaMemcpyHostToDevice );

	dim3 numThreadsPerBlock, numBlocks;
	numThreadsPerBlock.x = width;
//	numThreadsPerBlock.y 	= height;
	numBlocks.x = height;
//	numBlocks.y				= 1;

	cudaEvent_t start, end;
	cudaEventCreate( &start );
	cudaEventCreate( &end );
	float time;

	cout << "before gpu call" << endl;
	cudaEventRecord( start, 0 );
	calcLBPKernel<<< numBlocks, numThreadsPerBlock >>>( d_Src, d_Dst, width, height, mapping );
	cudaEventRecord( end, 0 );
	cudaEventSynchronize( end );
	cudaEventElapsedTime( &time, start, end );
	cudaEventDestroy( start );
	cudaEventDestroy( end );
	cout << "after gpu sync. Took " << time / 1000 << "s" << endl;

	cudaMemcpy( h_dst, d_Dst, sizeof(char) * height * width, cudaMemcpyDeviceToHost );
	cudaFree( d_Src );
	cudaFree( d_Dst );

}

bool cudaAvailable( void ) {
	int cnt;
	cudaGetDeviceCount( &cnt );
	cout << "GPU count: " << cnt << endl;
	return cnt > 0;
}

int main( int argc, char ** argv ) {

	if( !cudaAvailable() ) {
		return -1;
	}

	clock_t startTime, endTime;


#if 0
	unsigned char pixels[] = {78, 87, 84, 81, 92, 98,
				 75, 86, 82, 74, 82, 90,
				 77, 87, 85, 76, 74, 80,
				 91, 98, 91, 81, 77, 79,
				 90, 95, 85, 80, 84, 88,
				 91, 91, 83, 79, 86, 90};
	int w = 6, h = 6;
	cimg_library::CImg<unsigned char> image( pixels, w, h );

#else
	cimg_library::CImg<unsigned char> image( argv[1] );
	unsigned char * pixels = image.data();
	int w = image.width();
	int h = image.height();
#endif
	cout << "Image: " << image.pixel_type() << " (" << image.depth() << ") " << image.width() << "x"
				<< image.height() << ", # Channels=" << image.spectrum() << endl;

//	image.display();
//	for( int j = 0; j < image.height(); ++j ) {
//		for( int i = 0; i < image.width(); ++i ) {
//			printf( "%3d ", (int) image.atXY(i,j) );
//		}
//		printf( "\n" );
//	}

	LBPMapping mapping;

	startTime = clock();
	calcLBPGPU( pixels, pixels, w, h, &mapping );
	endTime = clock();
	// cudaDeviceReset must be called before exiting in order for profiling and
	    // tracing tools such as Nsight and Visual Profiler to show complete traces.
	cudaError_t cudaStatus = cudaDeviceReset();
	    if (cudaStatus != cudaSuccess) {
	        fprintf(stderr, "cudaDeviceReset failed!");
	        return 1;
	    }
#if 0
	for( int j = 0; j < h; j++ ) {
		for( int i = 0; i < w; ++i ) {
			printf( "%3d ", (int) pixels[j * w + i] );
		}
		printf( "\n" );
	}
#elseif 0
	image.assign( pixels, w, h );
		image.display();
#endif

	cout << "Example took " << double( endTime - startTime ) / double( CLOCKS_PER_SEC ) << "s"
				<< endl;



	return 0;

}
