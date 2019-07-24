#ifndef CGRID_H
#define CGRID_H

#include <iostream>
#include <cmath>
#include "Array2.h"
#include "Util.h"

enum CellType
{
	AIRCELL,FLUIDCELL, SOLIDCELL
};

class Grid {
public:
	Grid();
	Grid(float gravity_, int cell_nx, int cell_ny, float lx_);
	~Grid();

	unsigned int grid_width = 100;
	unsigned int grid_height = 100;

	//cell_nx, cell_ny: grid 전체 screen size
	// lx, ly : 한 cell 의 x,y 길이
	float gravity;
	float lx, ly;	//한 셀의 x,y 크기
	float h, overh;	// grid size a/n, 가로 셀 개수 n/a 
	
	//velocity
	Array2<float> u, v;
	Array2<float> du, dv; // differences for particle update
	Array2<CellType> marker; //identifies what sort of cell we have
	Array2<float> phi;	//distance field(phi)
	Array2<double> pressure;
	//
	Array2x3<float> poisson;
	Array2<double> preconditioner;
	Array2<double> m;
	Array2<double> r, z, s;

	void init(float gravity_, int cell_nx, int cell_ny, float lx_);
	float CFL();
	void bilerp_uv(float px, float py, float &pu, float &pv);
	void bary_x(float x, int &i, float &fx);
	void bary_y(float y, int &j, float &fy);
	void bary_x_centre(float x, int &i, float &fx);
	void bary_y_centre(float y, int &j, float &fy);
	void save_velocities();
	void add_gravity(float dt);
	void init_phi();
	void sweep_phi();
	void compute_distance_to_fluid();
	void extend_velocity();
	void sweep_velocity();
	void sweep_u(int i0, int i1, int j0, int j1);
	void sweep_v(int i0, int i1, int j0, int j1);
	void apply_boundary_conditions();
	void make_incompressible();
	void get_velocity_update();

	void find_divergence();
	void form_poisson();
	void form_preconditioner();
	void apply_poisson(const Array2<double>& x, Array2<double>& y);
	void apply_preconditioner(const Array2<double>& x, Array2<double>& y, Array2<double>& m);
	void solve_pressure(int maxits, double tolerance);
	void add_gradient();
	
private:
	
};
#endif