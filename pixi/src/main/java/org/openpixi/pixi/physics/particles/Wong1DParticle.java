package org.openpixi.pixi.physics.particles;

import org.openpixi.pixi.math.AlgebraElement;
import org.openpixi.pixi.math.ElementFactory;
import org.openpixi.pixi.math.GroupElement;
import org.openpixi.pixi.physics.RelativisticVelocity;

public class Wong1DParticle extends YangMillsParticle {

	public GroupElement U;
	public boolean updateCharge;
	public AlgebraElement E;

	public Wong1DParticle(int numberOfDimensions, int numberOfColors) {
		super(numberOfDimensions, numberOfColors);
		this.U = (new ElementFactory(numberOfColors)).groupIdentity();
		this.updateCharge = true;
	}

	public IParticle copy() {
		Wong1DParticle p = new Wong1DParticle(this.numberOfDimensions, this.numberOfColors);
		p.U = this.U.copy();

		for (int i = 0; i < this.numberOfDimensions; i++) {
			p.pos0[i] = this.pos0[i];
			p.pos1[i] = this.pos1[i];
			p.vel[i] = this.vel[i];
			p.acc[i] = this.acc[i];
		}

		p.setRadius(this.r);
		p.setDisplayColor(this.col);
		p.updateCharge = this.updateCharge;

		p.Q0 = this.Q0.copy();
		p.Q1 = this.Q1.copy();
		p.mass = this.mass;

		return p;
	}

	public double getGamma() {
		return Math.sqrt(1.0 + vel[0] * vel[0]);
	}

	public double getEnergy() {
		return mass * getGamma();
	}

}
