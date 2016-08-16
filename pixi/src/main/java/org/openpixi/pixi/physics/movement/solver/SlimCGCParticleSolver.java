package org.openpixi.pixi.physics.movement.solver;

import org.openpixi.pixi.math.AlgebraElement;
import org.openpixi.pixi.physics.force.Force;
import org.openpixi.pixi.physics.particles.IParticle;
import org.openpixi.pixi.physics.particles.SlimCGCParticle;

/**
 * Particle solver for slim CGC simulations. Basically the same as CGCParticleSolver.
 */
public class SlimCGCParticleSolver implements ParticleSolver {

	public void updatePosition(IParticle p, Force f, double dt) {
		SlimCGCParticle P = (SlimCGCParticle) p;

		for (int i = 0; i < P.getNumberOfDimensions(); i++) {
			P.setPosition(i, P.getPrevPosition(i) + P.getVelocity(i) * dt);
		}
	}


	public void updateCharge(IParticle p, Force f, double dt) {
		SlimCGCParticle P = (SlimCGCParticle) p;
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
