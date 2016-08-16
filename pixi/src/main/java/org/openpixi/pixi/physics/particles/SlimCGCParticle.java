package org.openpixi.pixi.physics.particles;

import org.openpixi.pixi.math.AlgebraElement;
import org.openpixi.pixi.math.ElementFactory;
import org.openpixi.pixi.math.GroupElement;

import java.awt.Color;
import java.io.Serializable;
import java.util.Arrays;

/**
 * Particle class to be used in CPIC simulations.
 */
public class SlimCGCParticle implements IParticle, Serializable {

	public AlgebraElement Q0;
	public AlgebraElement Q1;

    private double values[];
    private int dim;

    public int direction;
    public GroupElement U;
    public boolean updateCharge;

    // CONSTRUCTOR
	public SlimCGCParticle(int numberOfDimensions, int numberOfColors, int direction) {
		this.setNumberOfDimensions(numberOfDimensions);
		ElementFactory factory = new ElementFactory(numberOfColors);
		this.Q0 = factory.algebraZero();
		this.Q1 = factory.algebraZero();
        this.direction = direction;
        this.U = factory.groupIdentity();
        this.updateCharge = true;
	}

	// GETTERS

	public double getPosition(int i) {
		return values[dim + i];
	}

	public double getPrevPosition(int i) {
		return values[i];
	}

	public double getVelocity(int i) {
		return values[2*dim + i];
	}

	public double[] getPosition() {
	    return Arrays.copyOfRange(values, dim, 2*dim);
	}

	public double[] getPrevPosition() {
	    return Arrays.copyOfRange(values, 0, dim);
	}

	public double[] getVelocity() {
		return Arrays.copyOfRange(values, 2*dim, 3*dim);
	}

	public double getRadius() {
		return 0;
	}

	public Color getDisplayColor() {
		return Color.BLACK;
	}

	public int getNumberOfDimensions() {
		return this.dim;
	}

	// SETTERS

	public void setPosition(int i, double value) {
		this.values[dim + i] = value;
	}

	public void addPosition(int i, double value) {
		this.values[dim + i] += value;
	}

	public void setPrevPosition(int i, double value) {
		this.values[i] = value;
	}

	public void addPrevPosition(int i, double value) {
		this.values[i] += value;
	}

	public void setVelocity(int i, double value) {
		this.values[2*dim + i] = value;
	}

	public void addVelocity(int i, double value) {
		this.values[2*dim + i] += value;
	}

	public void setNumberOfDimensions(int numberOfDimensions) {
	    dim = numberOfDimensions;
	    // Store pos0, pos1, vel in one array.
        // index (0) - (DIM-1): pos0
        // index (DIM) - (2*DIM-1): pos1
        // index (2*DIM) - (3*DIM-1): vel
	    this.values = new double[3 * numberOfDimensions];
	}

	public void setRadius(double r) {

	}

	public void setDisplayColor(Color color) {

	}

	public void reassignValues() {
		double[] tempPos = Arrays.copyOfRange(values, 0, dim-1);
        System.arraycopy(values, dim, values, 0, dim);
        System.arraycopy(values, dim, tempPos, 0, dim);

		AlgebraElement tempQ = Q0;
		Q0 = Q1;
		Q1 = tempQ;
	}

	public IParticle copy() {
        int numberOfColors = 2;
		SlimCGCParticle p = new SlimCGCParticle(dim, numberOfColors, direction);

		for (int i = 0; i < dim; i++) {
		    p.setPosition(i, this.getPosition(i));
		    p.setPrevPosition(i, this.getPrevPosition(i));
		    p.setVelocity(i, this.getVelocity(i));
		}

		p.Q0 = this.Q0.copy();
		p.Q1 = this.Q1.copy();
        p.U = this.U.copy();
        p.updateCharge = this.updateCharge;

		return p;
	}

}
