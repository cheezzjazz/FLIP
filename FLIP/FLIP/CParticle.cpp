#include "CParticle.h"

Particle::Particle(Grid & grid_) : grid(grid_), num(0), sum(grid_.pressure.nx+1, grid_.pressure.ny+1)
{

}

Particle::~Particle()
{
}

void Particle::init(Grid & grid_)
{
	grid = grid_;
	num = 0;
	sum(grid_.pressure.nx + 1, grid_.pressure.ny + 1);
}

void Particle::add_particle(const Vec2<float>& px, const Vec2<float>& pu)
{
	x.push_back(px);
	u.push_back(pu);
	++num;
}

template<class T> 
void Particle::accumulate(T &accum, float q, int i, int j, float fx, float fy)
{
	float weight;

	weight = (1 - fx)*(1 - fy);
	accum(i, j) += weight*q;
	sum(i, j) += weight;

	weight = fx*(1-fy);
	accum(i+1, j) += weight*q;
	sum(i + 1, j) += weight;

	weight = (1 - fx)*fy;
	accum(i, j+1) += weight*q;
	sum(i, j+1) += weight;

	weight = fx*fy;
	accum(i+1, j+1) += weight*q;
	sum(i+1, j+1) += weight;
}

void Particle::move_particles_in_grid(float dt)//0.2*dt
{
	//a particle moves less than one grid cell in each substep
	Vec2<float> midx, gridu;
	float xmin = 1.001*grid.h, xmax = grid.lx - 1.001*grid.h;
	float ymin = 1.001*grid.h, ymax = grid.ly - 1.001*grid.h;
	for (int p = 0; p < num; p++)
	{
		grid.bilerp_uv(x[p][0], x[p][1], gridu[0], gridu[1]);
		midx = x[p] + 0.5*dt*gridu;
		clamp(midx[0], xmin, xmax);
		clamp(midx[1], ymin, ymax);
		//
		grid.bilerp_uv(midx[0], midx[1], gridu[0], gridu[1]);
		x[p] += dt*gridu;
		clamp(x[p][0], xmin, xmax);
		clamp(x[p][1], ymin, ymax);
	}
}

void Particle::transfer_to_grid()
{
	int p, i, ui, j, vj;
	float fx, ufx, fy, vfy;

	grid.u.zero();
	sum.zero();
	for (p = 0; p < num; p++)
	{
		grid.bary_x(x[p][0], ui, ufx);
		grid.bary_y_centre(x[p][1], j, fy);
		accumulate(grid.u, u[p][0], ui, j, ufx, fy);
	}
	for (j = 0; j < grid.u.ny; j++)
	{
		for (i = 0; i < grid.u.nx; i++)
		{
			if (sum(i, j) != 0) grid.u(i, j) /= sum(i, j);
		}
	}
	grid.v.zero();
	sum.zero();
	for (p = 0; p < num; p++)
	{
		grid.bary_x_centre(x[p][0], i, fx);
		grid.bary_y(x[p][1], vj, vfy);
		accumulate(grid.v, u[p][1], i, vj, fx, vfy);
	}
	for (j = 0; j < grid.v.ny; j++) for (i = 0; i < grid.v.nx; i++)
	{
		if (sum(i, j) != 0) grid.v(i, j) /= sum(i, j);
	}

	//identify where particles are in grid
	grid.marker.zero();
	for (p = 0; p < num; p++)
	{
		grid.bary_x(x[p][0], i, fx);
		grid.bary_y(x[p][1], j, fy);
		grid.marker(i, j) = FLUIDCELL;
	}
}

void Particle::update_from_grid()
{
	int p;
	int i, ui, j, vj;
	float fx, ufx, fy, vfy;
	for (p = 0; p < num; p++)
	{
		grid.bary_x(x[p][0], ui, ufx);
		grid.bary_x_centre(x[p][0], i, fx);
		grid.bary_y(x[p][1], vj, vfy);
		grid.bary_y_centre(x[p][1], j, fy);
		u[p] += Vec2<float>(grid.du.bilerp(ui, j, ufx, fy), grid.dv.bilerp(i, vj, fx, vfy)); //FLIP
	}
}