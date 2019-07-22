#ifndef PARTICLE_H
#define PARTICLE_H

#include "CGrid.h"
#include "Vec2.h"
#include "Util.h"
#include <iostream>
#include <vector>

class Particle {
public:
	Particle(Grid &grid_);
	~Particle();

	void init(Grid &grid_);
	Grid &grid;
	int num; // number of particles
	std::vector<Vec2<float>> x, u; //positions, and velocities

	//for transfer
	Array2<float> sum;

	void add_particle(const Vec2<float> &px, const Vec2<float> &pu);
	void move_particles_in_grid(float dt);
	void transfer_to_grid();

	void update_from_grid();

private:
	template<class T> void accumulate(T &accum, float q, int i, int j, float fx, float fy);
};
#endif