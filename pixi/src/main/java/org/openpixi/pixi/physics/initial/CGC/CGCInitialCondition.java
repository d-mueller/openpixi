package org.openpixi.pixi.physics.initial.CGC;

import org.openpixi.pixi.physics.Simulation;
import org.openpixi.pixi.physics.SimulationType;
import org.openpixi.pixi.physics.initial.IInitialCondition;

/**
 * This abstract class implements the common steps for CGC initial conditions:
 *      1) Initialize/load charge density. This depends on the particular CGC model.
 *      2) Solve the Poisson equation and initialize the fields in temporal gauge.
 *      3) Spawn particles based on the Gauss constraint and refine charge distribution.
 */
public class CGCInitialCondition implements IInitialCondition {

	/**
	 * Model for initial charge density.
	 */
	protected IInitialChargeDensity initialChargeDensity;

	/**
	 * CGC Poisson solver
	 */
	protected ICGCPoissonSolver solver;

	/**
	 * Particle creation algorithm.
	 */
	protected IParticleCreator initialParticleCreator;

	/**
	 * Option whether to compute the tadpole expectation value of the Wilson line.
	 */
	public boolean computeTadpole = false;

	/**
	 * File name or path for saving the tadpole output.
	 */
	public String tadpoleFilename = "tadpole.txt";

	/**
	 * Option whether to compute the dipole correlation function of the Wilson line.
	 */
	public boolean computeDipole = false;

	/**
	 * File name or path for saving the dipole output.
	 */
	public String dipoleFilename = "dipole.txt";


	public void initialize(Simulation s) {
		if(s.getSimulationType() == SimulationType.TemporalOptimizedCGCNGP) {
			// Reset all static objects needed for optimizations.
		}
	}

	/**
	 * Applies CGC initial conditions.
	 * @param s Reference to the Simulation object
	 */
	public void applyInitialCondition(Simulation s) {
		// Initialize charge density.
		initialChargeDensity.initialize(s);
		int direction = initialChargeDensity.getDirection();
		int orientation = initialChargeDensity.getOrientation();

		// Solve Poisson equation and set fields on the grid. Also computes Gauss constraint and saves it.
		solver.initialize(s);
		solver.solve(initialChargeDensity);

		// Compute tadpole, dipole, ...
		WilsonLineObservables wilson = new WilsonLineObservables(s, solver, initialChargeDensity);
		if(computeTadpole) {
			wilson.computeTadpole(tadpoleFilename);
		}
		if(computeDipole) {
			wilson.computeDipole(dipoleFilename);
		}

		// Spawn particles.
		switch (s.getSimulationType()) {
			case TemporalCGCNGP:
				initialParticleCreator = new LightConeNGPParticleCreator();
				break;
			case TemporalOptimizedCGCNGP:
				initialParticleCreator = new SuperCGCParticleCreator();
				break;
			default:
				System.out.println("CGCInitialCondition: simulation type not supported!");
				break;
		}

		initialParticleCreator.setGaussConstraint(solver.getGaussViolation());
		initialParticleCreator.initialize(s, direction, orientation);

		// Clear some memory.
		initialChargeDensity.clear();
	}

	public void setInitialChargeDensity(IInitialChargeDensity rho) {
		this.initialChargeDensity = rho;
	}

	public void setPoissonSolver(ICGCPoissonSolver solver) { this.solver = solver; }
}
