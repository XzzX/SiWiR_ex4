#ifndef VTK_HPP_INCLUDED
#define	VTK_HPP_INCLUDED

#include	<string>
#include	<iostream>
#include	<fstream>

#include	"format.h"

#include	"grid.hpp"
#include	"util.hpp"

using namespace std;

class	VTK{
private:
	const Params&	mParams;
	ofstream		mPVDStream;
public:
	const string PVD_FILENAME = "data/solution.pvd";
	const string VTU_FILENAME = "data/time_step_{:d}_{:d}.vtu";
	const string PVTU_FILENAME = "data/time_step_{:d}.pvtu";
	const string PVD_ENTRY = "  <DataSet timestep=\"{0:d}\" file=\"time_step_{0:d}.pvtu\"/>";
	
	VTK(const Params& params);
	~VTK();
	
	void addTimeStep(const int timestep, const Grid& u);
		
};

VTK::VTK(const Params& params)
	: mParams(params){
	if (mParams.rank == 0){
		mPVDStream.open(fmt::format(PVD_FILENAME));
		mPVDStream << "<?xml version=\"1.0\"?>" << endl;
		mPVDStream << "<VTKFile type=\"Collection\" version=\"0.1\">" << endl;
		mPVDStream << " <Collection>" << endl;
	}
}

VTK::~VTK(){
	if (mParams.rank == 0){
		mPVDStream << " </Collection>" << endl;
		mPVDStream << "</VTKFile>" << endl;
		mPVDStream.close();
	}
}

void VTK::addTimeStep(const int timestep, const Grid& u){
	if (mParams.rank == 0){
		//std::cout << timestep << std::endl;
		
		mPVDStream << fmt::format(PVD_ENTRY, timestep) << endl;
		
		ofstream fOut(fmt::format(PVTU_FILENAME, timestep));
		fOut << "<?xml version=\"1.0\"?>" << endl;
		fOut << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\">" << endl;
		fOut << " <PUnstructuredGrid>" << endl;
    
		fOut << "  <PPoints>" << endl;
		fOut << "   <PDataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\"/>" << endl;
		fOut << "  </PPoints>" << endl;
  
		fOut << "  <PPointData>" << endl;
		fOut << "   <PDataArray type=\"Float64\" Name=\"Temperature\" NumberOfComponents=\"1\" format=\"ascii\"/>" << endl;
		fOut << "  </PPointData>" << endl;
  
		for (int i = 0; i < mParams.size; ++i){
			fOut << fmt::format("  <Piece Source=\"time_step_{:d}_{:d}.vtu\"/>", timestep, i) << endl;
		}
  
		fOut << " </PUnstructuredGrid>" << endl;
		fOut << "</VTKFile>" << endl;
		fOut.close();
	}
	
	ofstream fOut(fmt::format(VTU_FILENAME, timestep, mParams.rank));
	fOut << "<?xml version=\"1.0\"?>" << endl;
	fOut << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">" << endl;
	fOut << " <UnstructuredGrid>" << endl;
	fOut << "  <Piece NumberOfPoints=\"" << mParams.bx * mParams.by << "\" NumberOfCells=\"" << mParams.bx * mParams.by <<"\">" << endl;
     
	fOut << "   <Points>" << endl;
	fOut << "    <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
	for (int y = 0; y < mParams.by; ++y){
		for (int x = 0; x < mParams.bx; ++x){
			fOut << "     " << mParams.getXCoord(x, y) << " " << mParams.getYCoord(x, y) << " 0.00" << std::endl;
		}
	}
	//	fOut << "     0.00 0.00 0.00" << endl;
	fOut << "    </DataArray>  " << endl;
	fOut << "   </Points>" << endl;
   
	fOut << "   <Cells>" << endl;
	fOut << "    <DataArray type=\"UInt32\" Name=\"connectivity\" format=\"ascii\">" << endl;
	for (int i = 0; i < mParams.bx * mParams.by; ++i){
		fOut << i << endl;
	}
	//fOut << "     0 1 2 3 4 5 6 7 8" << endl;
	fOut << "    </DataArray>" << endl;
	fOut << "    <DataArray type=\"UInt32\" Name=\"offsets\" format=\"ascii\">" << endl;
	for (int i = 1; i < mParams.bx * mParams.by + 1; ++i){
		fOut << i << endl;
	}
	//fOut << "     1 2 3 4 5 6 7 8 9" << endl;
	fOut << "    </DataArray>" << endl;
	fOut << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
	for (int i = 0; i < mParams.bx * mParams.by; ++i){
		fOut << 1 << endl;
	}
	//fOut << "     1 1 1 1 1 1 1 1 1" << endl;
	fOut << "    </DataArray>" << endl;
	fOut << "   </Cells>" << endl;
   
	fOut << "   <PointData>" << endl;
	fOut << "    <DataArray type=\"Float64\" Name=\"Temperature\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
	for (int y = 0; y < mParams.by; ++y){
		for (int x = 0; x < mParams.bx; ++x){
			fOut << "     "<< u(x, y) << std::endl;
		}
	}
	//fOut << "     0.0000 1.0000 2.0000" << endl;
	fOut << "    </DataArray>" << endl;
	fOut << "   </PointData>" << endl;
   
	fOut << "  </Piece>" << endl;
	fOut << " </UnstructuredGrid>" << endl;
	fOut << "</VTKFile>" << endl;

	fOut.close();
}

#endif
