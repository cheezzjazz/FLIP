#include <glad\glad.h>
#include <GLFW\glfw3.h>
#include <iostream>
#include "CParticle.h"
#include "CGrid.h"

class Solver {
public:
	Solver();
	Solver(int gridsize);
	~Solver();
		
	Grid *m_grid;
	Particle *m_particles;

	void init(int gridsize_);
	void init(int gridsize_, Grid & grid_, Particle & particles_);

	void update();
	void update(float * ptr, float * vertices);
	void stepPIC();
	void stepFLIP();
	void stepPICFLIP();

	float fluidphi(Grid &grid, float x, float y);
	void project(Grid & grid, float & x, float & y, float current, float target);
	void init_water_drop(Grid &grid, Particle &particles, int na, int nb);
	void advance_one_frame(Grid & grid, Particle & particles);
	void advance_one_frame(Grid &grid, Particle &particles, double frametime);
	void advanced_one_step(Grid &grid, Particle &particle, double dt);
private:
	unsigned int m_particle_count;
	//Particle *particles;
	//Grid *grid;
};