#ifndef GRID_HPP_INCLUDED
#define	GRID_HPP_INCLUDED

#include	<malloc.h>
#include	<iostream>
#include	<cstring>

/*
Matrix access object. Not responsible for memory management. Just handles access. Possibility for easy block matrix creation.
*/
class	Grid{
public:
	static const	int ALIGNMENT = 0x10;
	static const	int	PADDING = 10;

	///raw data for ROW MAJOR matrix
	double*		data;
	double*		backup;
	
	///number of rows
	int	sizeX;
	///number of cols
	int	sizeY;
	///size of border
	int borderSize;
	
	///size of leading dimension
	int	ld;

	///basic constructor
	Grid(const int _sizeX, const int _sizeY, const int _border = 1) : 
		sizeX(_sizeX + 2*_border), sizeY(_sizeY + 2*_border), borderSize(_border){
		
		ld = (sizeX + PADDING + 1) & 0xFFFFFFFE;
		backup = (double*) memalign(ALIGNMENT, sizeof(double) * ld * sizeY);
		data = backup + ld*borderSize + borderSize;
		//std::cout << ld << "\t" << borderSize << "\t" << backup << "\t" << data << std::endl;
		memset(backup, 0, sizeof(double) * ld * sizeY);
	}
	
	~Grid(){
		free(backup);
	}

	inline
	double& operator()(const int x, const int y) {
		return data[y * ld + x];
	}
	
	inline
	const double& operator()(const int x, const int y) const {
		return data[y * ld + x];
	}
};

#endif
