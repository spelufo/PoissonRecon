#include "PreProcessor.h"
#include "Reconstructors.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <random>
#include <iostream>
#include <fstream>
#include "MyMiscellany.h"

typedef unsigned int uint;
typedef Point<float, 3> Vec3f;

namespace PoissonVesuvius {

typedef float Real;
constexpr uint Dim = 3;
constexpr uint FEMSig = FEMDegreeAndBType<1, BOUNDARY_NEUMANN>::Signature;
typedef Reconstructor::Poisson::Implicit<Real, Dim, FEMSig> Implicit;

struct InputSampleStream : public Reconstructor::InputSampleStream<Real, Dim> {
	InputSampleStream(uint n, Vec3f *ps, Vec3f *ns) {
		_n = n; _ps = ps; _ns = ns;
		_i = 0;
	}
	void reset(void) { _i = 0; }
	bool base_read(Point<float, 3> &p, Point<float, 3> &n) {
		if (_i < _n) {
			p = _ps[_i];
			n = _ns[_i];
			// std::cout << "i: " << _i << std::endl;
			// std::cout << "p: " << p << std::endl;
			// std::cout << "n: " << n << std::endl;
			_i++;
			return true;
		}
		return false;
	}
protected:
	uint _n; Vec3f *_ps; Vec3f *_ns;
	uint _i;
};








// A stream for generating random samples on the sphere
struct SphereSampleStream : public Reconstructor::InputSampleStream< Real , Dim >
{
	// from https://en.cppreference.com/w/cpp/numeric/random/uniform_real_distribution
	std::random_device randomDevice;
	std::default_random_engine generator;
	std::uniform_real_distribution< Real > distribution;

	// Constructs a stream that contains the specified number of samples
	SphereSampleStream( unsigned int sz ) : _size(sz) , _current(0) , generator(0) , distribution((Real)-1.0,(Real)1.0) {}

	// Overrides the pure abstract method from InputSampleStream< Real , Dim >
	void reset( void ){ generator.seed(0) ; _current = 0; }

	// Overrides the pure abstract method from InputSampleStream< Real , Dim >
	bool base_read( Point< Real , Dim > &p , Point< Real , Dim > &n )
	{
		if( _current<_size )
		{
			p = n = RandomSpherePoint( generator , distribution );
			_current++;
			return true;
		}
		else return false;
	}

	static Point< Real , Dim > RandomSpherePoint( std::default_random_engine &generator , std::uniform_real_distribution< Real > &distribution )
	{
		while( true )
		{
			Point< Real , Dim > p;
			for( unsigned int d=0 ; d<Dim ; d++ ) p[d] = distribution( generator );
			if( Point< Real , Dim >::SquareNorm( p )<1 ) return p / (Real)sqrt( Point< Real , Dim >::SquareNorm(p) );
		}
	}
protected:
	unsigned int _size , _current;
};





struct PolygonStream : public Reconstructor::OutputFaceStream<2> {
	PolygonStream(std::vector<std::vector<int>> &polygonStream) : _polygons(polygonStream) {}

	void base_write(const std::vector<node_index_type> &polygon) {
		std::vector<int> poly(polygon.size());
		for (uint i=0; i < polygon.size(); i++) poly[i] = (int)polygon[i];
		_polygons.push_back(poly);
	}

protected:
	std::vector<std::vector<int>> &_polygons;
};

struct VertexStream : public Reconstructor::OutputVertexStream<Real, Dim> {
	VertexStream(std::vector<Real> &vCoordinates) : _vCoordinates(vCoordinates) {}

	void base_write(Point<Real, Dim> p, Point<Real, Dim>, Real){
		for (uint d=0; d < Dim; d++) _vCoordinates.push_back(p[d]);
	}
protected:
	std::vector<Real> &_vCoordinates;
};

// ASCII :(
void WritePly(std::string fileName, size_t vNum, const Real *vCoordinates, const std::vector<std::vector<int>> &polygons) {
	std::fstream file(fileName, std::ios::out);
	file << "ply" << std::endl;
	file << "format ascii 1.0" << std::endl;
	file << "element vertex " << vNum << std::endl;
	file << "property float x" << std::endl << "property float y" << std::endl << "property float z" << std::endl;
	file << "element face " << polygons.size() << std::endl;
	file << "property list uchar int vertex_indices" << std::endl;
	file << "end_header" << std::endl;

	for(size_t i=0 ; i<vNum ; i++) {
		file << vCoordinates[3*i+0] << " " << vCoordinates[3*i+1] << " " << vCoordinates[3*i+2];
		file << std::endl;
	}
	for(const auto &polygon : polygons) {
		file << polygon.size();
		for(auto vIdx : polygon) file << " " << vIdx;
		file << std::endl;
	}
}

void PoissonRecon(
	uint n, Vec3f *ps, Vec3f *ns,
	const char* outPath,
	Reconstructor::Poisson::SolutionParameters<Real> solverParams,
	Reconstructor::LevelSetExtractionParameters extractionParams
) {
	#ifdef _OPENMP
		ThreadPool::Init(ThreadPool::OPEN_MP, std::thread::hardware_concurrency());
	#else
		ThreadPool::Init(ThreadPool::THREAD_POOL, std::thread::hardware_concurrency());
	#endif

	InputSampleStream sampleStream(n, ps, ns);
	// SphereSampleStream sampleStream(1000);
	// Reconstructor::Poisson::EnvelopeMesh<Real, Dim> *envelopeMesh = NULL; // TODO
	Implicit implicit(sampleStream, solverParams);

	// TODO: Output in memory.
	std::vector<std::vector<int>> polygons;
	std::vector<Real> vCoordinates;
	VertexStream vStream(vCoordinates);
	PolygonStream pStream(polygons);
	implicit.extractLevelSet(vStream, pStream, extractionParams);

	WritePly(outPath, vStream.size(), vCoordinates.data(), polygons);

	ThreadPool::Terminate();
}

} // PoissonVesuvius

extern "C"
void poisson_recon(
	uint n, Vec3f *ps, Vec3f *ns,
	const char* outPath,
	Reconstructor::Poisson::SolutionParameters<float> solverParams,
	Reconstructor::LevelSetExtractionParameters extractionParams
) {
	PoissonVesuvius::PoissonRecon(n, ps, ns, outPath, solverParams, extractionParams);
}
