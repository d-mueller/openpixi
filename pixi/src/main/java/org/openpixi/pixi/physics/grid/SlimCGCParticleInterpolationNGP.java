package org.openpixi.pixi.physics.grid;

import org.openpixi.pixi.math.AlgebraElement;
import org.openpixi.pixi.math.GroupElement;
import org.openpixi.pixi.physics.particles.CGCParticle;
import org.openpixi.pixi.physics.particles.IParticle;
import org.openpixi.pixi.physics.particles.SlimCGCParticle;
import org.openpixi.pixi.physics.util.GridFunctions;

/**
 * This interpolation algorithm is used for CGC simulations in the lab frame. The particles in this type of simulation
 * merely act as 'static' sources for the current. It is assumed that this kind of particle moves along a grid axis such
 * that there is no ambiguity in defining parallel transport for the color charges of the particles.
 */
public class SlimCGCParticleInterpolationNGP implements  InterpolatorAlgorithm {
	public void interpolateToGrid(IParticle p, Grid g) {
		SlimCGCParticle P = (SlimCGCParticle) p;

		double at = g.getTemporalSpacing();
		double as = g.getLatticeSpacing();
		int direction = P.direction;

		// Particle positions
		double[] oldPosition = P.getPrevPosition();
		double[] newPosition = P.getPosition();

		// check if one cell or two cell move
		int[] ngpOld = GridFunctions.nearestGridPoint(oldPosition, as);
		int[] ngpNew = GridFunctions.nearestGridPoint(newPosition, as);

		if(ngpOld[direction] == ngpNew[direction]) {
			// one cell move
		} else {
			// two cell move
			int cellIndexOld = g.getCellIndex(ngpOld);
			int cellIndexNew = g.getCellIndex(ngpNew);

			if(P.getVelocity(direction) > 0) {
				AlgebraElement J = P.Q0.mult(as / at);
				g.addJ(cellIndexOld, direction, J);
			} else {
				AlgebraElement J = P.Q1.mult(- as / at);
				g.addJ(cellIndexNew, direction, J);
			}
		}
	}

	public void interpolateChargedensity(IParticle p, Grid g) {
		SlimCGCParticle P = (SlimCGCParticle) p;

		double as = g.getLatticeSpacing();

		// "Floored" grid points of the particle
		int[] gridPosOld = GridFunctions.nearestGridPoint(P.getPrevPosition(), as);

		// Cell indices
		int cellIndexOld = g.getCellIndex(gridPosOld);

		g.addRho(cellIndexOld, P.Q0);

	}

	public void interpolateToParticle(IParticle p, Grid g) {
		// Compute parallel transport for the particle.
		SlimCGCParticle P = (SlimCGCParticle) p;
		double as = g.getLatticeSpacing();
		int direction = P.direction;

		// Particle positions
		double[] oldPosition = P.getPrevPosition();
		double[] newPosition = P.getPosition();

		// check if one cell or two cell move
		int[] ngpOld = GridFunctions.nearestGridPoint(oldPosition, as);
		int[] ngpNew = GridFunctions.nearestGridPoint(newPosition, as);

		if(ngpOld[direction] == ngpNew[direction]) {
			// one cell move: do nothing!
			//P.U = identity;
		} else {
			// two cell move
			int cellIndexOld = g.getCellIndex(ngpOld);
			int cellIndexNew = g.getCellIndex(ngpNew);

			if(P.getVelocity(direction) > 0) {
				P.U = g.getUnext(cellIndexOld, direction);
			} else {
				P.U = g.getUnext(cellIndexNew, direction).adj();
			}
			P.updateCharge = true;
		}
	}
}
