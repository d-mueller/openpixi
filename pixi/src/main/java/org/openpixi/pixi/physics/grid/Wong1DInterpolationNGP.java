package org.openpixi.pixi.physics.grid;

import org.openpixi.pixi.math.AlgebraElement;
import org.openpixi.pixi.physics.particles.IParticle;
import org.openpixi.pixi.physics.particles.Wong1DParticle;
import org.openpixi.pixi.physics.util.GridFunctions;

public class Wong1DInterpolationNGP implements  InterpolatorAlgorithm {
	public void interpolateToGrid(IParticle p, Grid g) {
		Wong1DParticle P = (Wong1DParticle) p;

		double at = g.getTemporalSpacing();
		double as = g.getLatticeSpacing();

		// Particle positions
		double[] oldPosition = P.pos0;
		double[] newPosition = P.pos1;

		// check if one cell or two cell move
		int[] ngpOld = GridFunctions.nearestGridPoint(oldPosition, as);
		int[] ngpNew = GridFunctions.nearestGridPoint(newPosition, as);

		if(ngpOld[0] == ngpNew[0]) {
			// one cell move
		} else {
			// two cell move
			int cellIndexOld = g.getCellIndex(ngpOld);
			int cellIndexNew = g.getCellIndex(ngpNew);

			if(P.vel[0] > 0) {
				AlgebraElement J = P.Q0.mult(as / at);
				g.addJ(cellIndexOld, 0, J);
			} else {
				AlgebraElement J = P.Q1.mult(- as / at);
				g.addJ(cellIndexNew, 0, J);
			}
		}
	}

	public void interpolateChargedensity(IParticle p, Grid g) {
		Wong1DParticle P = (Wong1DParticle) p;

		double as = g.getLatticeSpacing();

		// "Floored" grid points of the particle
		int[] gridPosOld = GridFunctions.nearestGridPoint(P.pos0, as);

		// Cell indices
		int cellIndexOld = g.getCellIndex(gridPosOld);

		g.addRho(cellIndexOld, P.Q0);

	}

	public void interpolateToParticle(IParticle p, Grid g) {
		// Compute parallel transport for the particle.
		Wong1DParticle P = (Wong1DParticle) p;
		double as = g.getLatticeSpacing();

		// Particle positions
		double[] oldPosition = P.pos0;
		double[] newPosition = P.pos1;

		// check if one cell or two cell move
		int[] ngpOld = GridFunctions.nearestGridPoint(oldPosition, as);
		int[] ngpNew = GridFunctions.nearestGridPoint(newPosition, as);

		int cellIndexOld = g.getCellIndex(ngpOld);
		int cellIndexNew = g.getCellIndex(ngpNew);

		if(ngpOld[0] == ngpNew[0]) {
			// one cell move: do nothing!
			//P.U = identity;
		} else {
			// two cell move

			if(P.vel[0] > 0) {
				P.U = g.getUnext(cellIndexOld, 0);
			} else {
				P.U = g.getUnext(cellIndexNew, 0).adj();
			}
			P.updateCharge = true;
		}

		// Interpolate chromoelectric field
		P.E = g.getE(cellIndexNew, 0);
	}
}
