#include "CSolver.h"

Solver::Solver()
{
}

Solver::Solver(int gridsize)
{
	init(gridsize);
}

Solver::~Solver()
{
	if(m_particles != NULL)
		delete m_particles;
	if(m_grid != NULL)
		delete m_grid;
}

void Solver::init(int gridsize_)
{
	m_grid = new Grid(9.8, gridsize_, gridsize_, 1);
	m_particles = new Particle(*m_grid);

	init_water_drop(*m_grid, *m_particles, 2, 2);
}

void Solver::init(int gridsize_, Grid &grid_, Particle &particles_)
{
	Grid grid(9.8, gridsize_, gridsize_, 1);
	Particle particles(grid);
	init_water_drop(grid, particles, 2, 2);
}

void Solver::update()
{
}

void Solver::update(float *ptr, float * vertices)
{
	if (!ptr || !vertices)
		return;

	stepFLIP();

	for (int p = 0; p < m_particle_count; p++)
	{
		*vertices = m_particles->x[p][0];
		*ptr = *vertices;
		std::cout << "[x," << p << "] = " << *ptr << "," << *vertices << std::endl; ++ptr; ++vertices;
		*vertices = m_particles->x[p][1];
		*ptr = *vertices;
		std::cout << "[y," << p << "] = " << *ptr << "," << *vertices << std::endl; ++ptr; ++vertices;
	}
}

void Solver::stepPIC()
{
}

void Solver::stepFLIP()
{
	advance_one_frame(*m_grid, *m_particles, 1./30);
}

void Solver::stepPICFLIP()
{
}

float Solver::fluidphi(Grid & grid, float x, float y)
{
	return min(y-0.05*grid.ly, sqrt(sqr(x-0.5*grid.lx)+sqr(y-0.5*grid.ly))-0.05*grid.lx);
}

void Solver::project(Grid &grid, float &x, float &y, float current, float target)
{
	float dpdx = (fluidphi(grid, x+1e-4, y) - fluidphi(grid, x-1e-4, y)) / 2e-4;
	float dpdy = (fluidphi(grid, x, y+1e-4) - fluidphi(grid, x, y-1e-4)) / 2e-4;
	float scale = (target - current) / sqrt(dpdx*dpdx + dpdy*dpdy);
	x += scale*dpdx;
	y += scale*dpdy;
}

void Solver::init_water_drop(Grid & grid, Particle & particles, int na, int nb)
{
	float x, y, phi;
	
	for (int i = 1; i < grid.marker.nx - 1; i++)	//except for boundary
	{
		for (int j = 1; j < grid.marker.ny - 1; j++)
		{
			for (int a = 0; a < na; ++a)
			{
				for (int b = 0; b < nb; ++b)
				{
					x = (i + (a + 0.1f + 0.8f*rand() / (float)RAND_MAX) / na) * grid.h; //0.0~0.8 random, random seeding
					y = (j + (b + 0.1f + 0.8f*rand() / (float)RAND_MAX) / na) * grid.h;
					phi = fluidphi(grid, x, y);
					if (phi > -0.25*grid.h / na)
						continue;
					else if (phi > -1.5 * grid.h / na)
					{
						project(grid, x, y, phi, -0.75*grid.h / na);
						phi = fluidphi(grid, x, y);
						project(grid, x, y, phi, -0.75*grid.h / na);
						phi = fluidphi(grid, x, y);
					}
					particles.add_particle(Vec2<float>(x, y), Vec2<float>(0, 0));
				}
			}
		}
	}
}

void Solver::advance_one_frame(Grid & grid, Particle & particles, double frametime)
{
	m_particle_count = 0;
	double t = 0;
	double dt;
	bool finished = false;
	while (!finished)
	{
		dt = 2 * grid.CFL();
		if (t + dt >= frametime)
		{
			dt = frametime - t;
			finished = true;
			m_particle_count = particles.num;
		}
		else if (t + 1.5*dt >= frametime)
		{
			dt = 0.5*(frametime - t);
			std::cout << "advancing " << dt << "(to " << 100 * (t + dt) / frametime << "of frame)" << std::endl;
			advanced_one_step(grid, particles, dt);
			t += dt;
		}
	}
}

void Solver::advanced_one_step(Grid & grid, Particle & particle, double dt)
{
	for (int i = 0; i < 5; i++)
		particle.move_particles_in_grid(0.2*dt);

	particle.transfer_to_grid();
	grid.save_velocities();
	grid.add_gravity(dt);
	grid.compute_distance_to_fluid();
	grid.extend_velocity();
	grid.apply_boundary_conditions();
	grid.make_incompressible();
	grid.extend_velocity();
	grid.get_velocity_update();
	particle.update_from_grid();
}
