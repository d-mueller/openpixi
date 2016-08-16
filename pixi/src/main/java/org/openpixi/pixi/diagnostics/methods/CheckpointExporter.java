package org.openpixi.pixi.diagnostics.methods;

import org.openpixi.pixi.diagnostics.Diagnostics;
import org.openpixi.pixi.physics.Simulation;
import org.openpixi.pixi.physics.grid.Grid;
import org.openpixi.pixi.physics.particles.IParticle;

import java.io.IOException;
import java.util.ArrayList;

public class CheckpointExporter implements Diagnostics {

	public void initialize(Simulation simulation) {

	}

	public void calculate(Grid grid, ArrayList<IParticle> particles, int steps) throws IOException {
		// 1) Write header with meta info:
		// number of dimensions, number of colors, grid size, spatial and temporal lattice spacing, coupling constant
		// number of particles

		/*
			Note: loading checkpoints is implemented as an initial condition, which means that the grid will already be
			initialized. The header is therefore is a bit meaningless, but we can use it to check for incompatibilities
			before loading the whole file.
		 */

		// 2) Write U and Unext data.

		// 3) Write E data.

		// 4) Write particle data (2 position vectors, velocity, 2 charge elements and the direction integer)

	}
}
