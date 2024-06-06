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

// Julia callbacks to build output mesh.
typedef void (*callback_vert_new)(float, float, float, float, float, float, float);
typedef void (*callback_tri_new)(uint32_t, uint32_t, uint32_t);
typedef void (*callback_poly_new)(uint64_t);
typedef void (*callback_poly_index)(uint32_t);

namespace PoissonVesuvius {

typedef float Real;
constexpr uint Dim = 3;
constexpr uint FEMSig = FEMDegreeAndBType<1, BOUNDARY_NEUMANN>::Signature;
typedef Reconstructor::Poisson::Implicit<Real, Dim, FEMSig> Implicit;

struct InputSampleStream : public Reconstructor::InputSampleStream<Real, Dim> {
	InputSampleStream(uint n, Vec3f *ps, Vec3f *ns) : _n(n), _ps(ps), _ns(ns), _i(0) {}
	void reset(void) { _i = 0; }
	bool base_read(Point<float, 3> &p, Point<float, 3> &n) {
		if (_i < _n) {
			p = _ps[_i]; n = _ns[_i]; _i++;
			return true;
		}
		return false;
	}
protected:
	uint _n; Vec3f *_ps; Vec3f *_ns; uint _i;
};

// Adapter streams to output to julia callbacks, so that julia does the alloc.
struct VertexStream : public Reconstructor::OutputVertexStream<Real, Dim> {
	VertexStream(callback_vert_new cbvert) : _cbvert(cbvert) {}
	void base_write(Point<Real, Dim> p, Point<Real, Dim> n, Real density){
		_cbvert(p[0], p[1], p[2], n[0], n[1], n[2], density);
	}
protected:
	callback_vert_new _cbvert;
};
struct PolygonStream : public Reconstructor::OutputFaceStream<2> {
	PolygonStream(callback_tri_new cbtri, callback_poly_new cbpoly, callback_poly_index cbpolyidx)
		: _cbtri(cbtri), _cbpoly(cbpoly), _cbpolyidx(cbpolyidx) {}
	void base_write(const std::vector<node_index_type> &polygon) {
		auto sides = polygon.size();
		if (sides == 3) {
			_cbtri(polygon[0], polygon[1], polygon[2]);
		} else {
			_cbpoly(sides);
			for (uint i=0; i < polygon.size(); i++) _cbpolyidx(polygon[i]);
		}
	}
protected:
	callback_tri_new    _cbtri;
	callback_poly_new   _cbpoly;
	callback_poly_index _cbpolyidx;
};

// // ASCII :(
// void WritePly(std::string fileName, size_t vNum, const Real *vCoordinates, const std::vector<std::vector<int>> &polygons) {
// 	std::fstream file(fileName, std::ios::out);
// 	file << "ply" << std::endl;
// 	file << "format ascii 1.0" << std::endl;
// 	file << "element vertex " << vNum << std::endl;
// 	file << "property float x" << std::endl << "property float y" << std::endl << "property float z" << std::endl;
// 	file << "element face " << polygons.size() << std::endl;
// 	file << "property list uchar int vertex_indices" << std::endl;
// 	file << "end_header" << std::endl;
//
// 	for(size_t i=0 ; i<vNum ; i++) {
// 		file << vCoordinates[3*i+0] << " " << vCoordinates[3*i+1] << " " << vCoordinates[3*i+2];
// 		file << std::endl;
// 	}
// 	for(const auto &polygon : polygons) {
// 		file << polygon.size();
// 		for(auto vIdx : polygon) file << " " << vIdx;
// 		file << std::endl;
// 	}
// }

void PoissonRecon(
	uint n, Vec3f *ps, Vec3f *ns,
	callback_vert_new cbvert, callback_tri_new cbtri, callback_poly_new cbpoly, callback_poly_index cbpolyidx,
	Reconstructor::Poisson::SolutionParameters<Real> solverParams,
	Reconstructor::LevelSetExtractionParameters extractionParams
) {

	// By default when I first compiled this it was using openmp. I started seeing
	// hangs when calling from julia, so I changed this not to use openmp. If the
	// problem is dlclose we can stop doing that too. Which thread pool is faster?
	// https://github.com/JuliaLang/julia/issues/10938#issuecomment-95001415
	ThreadPool::Init(ThreadPool::THREAD_POOL, std::thread::hardware_concurrency());

	InputSampleStream sampleStream(n, ps, ns);
	// Reconstructor::Poisson::EnvelopeMesh<Real, Dim> *envelopeMesh = NULL;
	Implicit implicit(sampleStream, solverParams);

	// TODO: Output in memory.
	VertexStream vStream(cbvert);
	PolygonStream pStream(cbtri, cbpoly, cbpolyidx);
	implicit.extractLevelSet(vStream, pStream, extractionParams);

	// WritePly(outPath, vStream.size(), vCoordinates.data(), polygons);

	ThreadPool::Terminate();
}

} // PoissonVesuvius

extern "C"
void poisson_recon(
	uint n, Vec3f *ps, Vec3f *ns,
	callback_vert_new cbvert, callback_tri_new cbtri, callback_poly_new cbpoly, callback_poly_index cbpolyidx,
	Reconstructor::Poisson::SolutionParameters<float> solverParams,
	Reconstructor::LevelSetExtractionParameters extractionParams
) {
	PoissonVesuvius::PoissonRecon(n, ps, ns, cbvert, cbtri, cbpoly, cbpolyidx, solverParams, extractionParams);
}
