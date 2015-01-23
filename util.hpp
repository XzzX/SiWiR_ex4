#ifndef	UTIL_HPP_INCLUDED
#define	UTIL_HPP_INCLUDED

#include	<stdlib.h>
#include	<iostream>
#include	<sstream>
#include	<fstream>
#include	<cmath>
#include	<limits>       // std::numeric_limits

/**
  Converts a string to an arbitrary type. >> operator must be defined for the target type.
  @param string string which should be converted
  @return converted string
  **/
template<typename T>
T StringTo(const std::string& string) {
	T valor;

	std::stringstream stream(string);
	stream >> valor;
	return valor;
}

constexpr	double f(const double x, const double y) {
	return 4.0*M_PI*M_PI*sin(2.0*M_PI*x)*sinh(2.0*M_PI*y);
}

double initial(const double x, const double y) {
	//return sin(M_PI*x)*sin(M_PI*y);
	double tempX = x - 0.5;
	double tempY = y - 0.5;
	return sin(2*M_PI* sqrt(tempX*tempX + tempY*tempY));
}

class	Params{
public:
	static constexpr double k  = 2.0 * M_PI;
	
	static const int UP        = 0;
	static const int DOWN      = 1;
	static const int LEFT      = 2;
	static const int RIGHT     = 3;
	static const int UPLEFT    = 4;
	static const int UPRIGHT   = 5;
	static const int DOWNLEFT  = 6;
	static const int DOWNRIGHT = 7;
	
	bool	output             = false;
	bool	debugOutput        = false;
	
	int nx                     = 0;
	int	ny                     = 0;
	int	c                      = 0;
	double eps                 = 0.0;
	double eps2                = 0.0;
	int timesteps				= 0;
	double tau					= 0.0;
	double kappa				= 0.0;
	double alpha				= 0.0;
	int vtk_spacing				= 0;

	double	hx                 = 0.0;
	double	hy                 = 0.0;

	double	invHx2             = 0.0;
	double	invHy2             = 0.0;

	double	preF               = 0.0;
	
	int	size                   = 0; // The total number of processes
	int rank                   = 0; // The rank/number of this process (within MPI_COMM_WORLD)
	
	// Definition and initialization of the dimension array 'dims'
	// In this case, we will construct a 2D-grid with 2 processes in x-
	// and 2 processes in y-direction
	int dims[2]                = {2,2};

	// Definition and initialization of the behavior of the Cartesian topology
	int periods[2]             = {0,0}; // not periodic!
	const int reorder          = 1;  // allow reordering of process ranks
	
	// The rank/number of this process (within the Cartesian topology)
	int cartrank               = 0;

	// Coordinates of this process within the Cartesian topology
	int coords[2]              = {0,0};

	// Definition of the neighbors array 'nrbs'
	// Will contain the ranks of the neighbors in the order UP, DOWN, LEFT and RIGHT
	int nbrs[8]                = {0,0,0,0,0,0,0,0};
	
	//size of block
	int	bx	                   = 0;
	int	by	                   = 0;
	int offsetX                = 0;
	int offsetY                = 0;
		
	Params(int argc, char **argv){
		if (argc < 10) {
			std::cout << "Invalid number of arguments!" << std::endl;
			std::cout << "./rgbs nx ny c eps" << std::endl;
			exit(EXIT_FAILURE);
		}
		nx          = StringTo<int>(argv[1]);
		ny          = StringTo<int>(argv[2]);
		c           = StringTo<int>(argv[3]);
		eps         = StringTo<double>(argv[4]);
		eps2 		= copysign(eps * eps, eps);
		timesteps	= StringTo<int>(argv[5]);;
		tau			= StringTo<double>(argv[6]);
		kappa		= StringTo<double>(argv[7]);
		alpha		= StringTo<double>(argv[8]);
		vtk_spacing = StringTo<int>(argv[9]);

		// output configuration parameters
		/*
		if (rank==0){
			std::cout << "nx," << nx << std::endl;
			std::cout << "ny," << ny << std::endl;
			std::cout << "c," << c << std::endl;
			std::cout << "eps," << eps << std::endl;
			std::cout << "eps2," << eps2 << std::endl;
		}
		*/
		
		// calculate global parameters
		hx     = 1.0 / nx;
		hy     = 1.0 / ny;

		invHx2 = 1.0 / hx / hx;
		invHy2 = 1.0 / hy / hy;

		preF   = ( 2 * invHx2 + 2 * invHy2 );
		
		//make real number of points
		++nx;
		++ny;
		
		for (int i = 5; i < argc; ++i){
			std::string cmd = argv[i];
			
			if (cmd.compare("-v") == 0)
				output = true;
				
			if (cmd.compare("-vv") == 0)
				debugOutput = true;
		}
	}
	
	void	subdivideGrid(){
		dims[0] = 1;
		dims[1] = 1;
		int	tempSize = size;
		while (((tempSize&0x1) == 0) && (tempSize!=0)){
			tempSize /= 2;
			if (tempSize > dims[0]) {
				dims[0] *= 2;
			}
		}
		dims[1] = size / dims[0];
		if (rank == 0){
			std::cout << "blocks," << dims[0] << "x" << dims[1] << std::endl;
		}
	}
	
	void	createBlock(){
		int	reminder = (nx - 2) % dims[0];
		bx      = (nx - 2) / dims[0];
		offsetX = bx * coords[0] + 1;
		if (coords[0] < reminder) {
			++bx;
			offsetX += coords[0];
		} else {
			offsetX += reminder;
		}
		
		reminder = (ny - 2) % dims[1];
		by      = (ny - 2) / dims[1];
		offsetY = by * coords[1] + 1;
		if (coords[1] < reminder) {
			++by;
			offsetY += coords[1];
		} else {
			offsetY += reminder;
		}
		
		for (int i = 0; i < size; ++i){
			if (i == rank)
				std::cout << "block," << rank << "\t" << coords[0] << "\t" << coords[1] << "\t" << offsetX << "\t" << offsetY << "\t" << bx << "\t" << by << std::endl;
			MPI_Barrier( MPI_COMM_WORLD );
		}
	}
	
	inline	
	int	getBottomBorderOffset(){
		return (coords[1] == 0) ? -1 : 0;
	}
	
	inline	
	int	getTopBorderOffset(){
		return (coords[1] == dims[1]-1) ? 1 : 0;
	}
	
	inline	
	int	getLeftBorderOffset(){
		return (coords[0] == 0) ? -1 : 0;
	}
	
	inline	
	int	getRightBorderOffset(){
		return (coords[0] == dims[0]-1) ? 1 : 0;
	}
	
	inline
	double	getXCoord(const int indexX, const int indexY) const{
		(void)(indexY);
		return	hx * (indexX + offsetX);
	}
	
	inline
	double	getYCoord(const int indexX, const int indexY) const{
		(void)(indexX);
		return	hy * (indexY + offsetY);
	}

};


	
#endif