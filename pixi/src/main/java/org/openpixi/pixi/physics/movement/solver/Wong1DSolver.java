package org.openpixi.pixi.physics.movement.solver;

import org.openpixi.pixi.math.AlgebraElement;
import org.openpixi.pixi.physics.force.Force;
import org.openpixi.pixi.physics.grid.Wong1DInterpolationNGP;
import org.openpixi.pixi.physics.particles.CGCParticle;
import org.openpixi.pixi.physics.particles.IParticle;
import org.openpixi.pixi.physics.particles.Wong1DParticle;

/**
 * Particle solver for CGC simulations. The particles stay on fixed trajectories moving at light speed. The charge has
 * to be parallel transporting along the trajectory. It is assumed that particles move on a grid line such that parallel
 * transport using gauge links is well defined.
 */
public class Wong1DSolver implements ParticleSolver {

	public void updatePosition(IParticle p, Force f, double dt) {
		Wong1DParticle P = (Wong1DParticle) p;
		P.pos1[0] = P.pos0[0] + P.vel[0] / P.getGamma() * dt;
	}

	public void updateVelocity(IParticle p, Force f, double step)
	{
		Wong1DParticle P = (Wong1DParticle) p;
		P.vel[0] += f.getForce(0, p) / P.mass * step;
	}

	public void updateCharge(IParticle p, Force f, double dt) {
		Wong1DParticle P = (Wong1DParticle) p;
		if(P.updateCharge) {
			// Charge has to be parallel transported.
			P.Q1 = P.Q0.act(P.U.adj());
			P.updateCharge = false;
		} else {
			// No update needed, just switch Q1 and Q0.
			AlgebraElement Q = P.Q1;
			P.Q1 = P.Q0;
			P.Q0 = Q;
		}
	}

	public void prepare(IParticle p, Force f, double step) {
		// Not implemented.
	}

	public void complete(IParticle p, Force f, double step) {
		// Not implemented.
	}
}
