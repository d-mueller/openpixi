package org.openpixi.pixi.physics.initial.YM;

import org.openpixi.pixi.math.AlgebraElement;
import org.openpixi.pixi.math.GroupElement;
import org.openpixi.pixi.math.SU2AlgebraElement;
import org.openpixi.pixi.math.SU2GroupElement;
import org.openpixi.pixi.physics.Simulation;
import org.openpixi.pixi.physics.initial.CGC.FourierFunctions;
import org.openpixi.pixi.physics.initial.IInitialCondition;
import org.openpixi.pixi.physics.util.GridFunctions;

import java.util.Random;

public class SU2GlasmaFluxTubes implements IInitialCondition {
	private int direction;
	private double mu;
	private double IR;
	private double UVT;
	private Random rand;
	private int totalCells;
	private int[] transNumCells;
	private int totalTransCells;
	private int numberOfColors;
	private int numberOfComponents;
	private double aT;
	private double g;
	private Simulation s;

	public SU2GlasmaFluxTubes(int direction, double mu, double IR, double UVT, boolean useSeed, int seed) {
		this.direction = direction;
		this.mu = mu;
		this.IR = IR;
		this.UVT = UVT;

		this.rand = new Random();
		if (useSeed) {
			rand.setSeed(seed);
		}
	}

	/**
	 * Initializes the Glasma fields according to the MV model.
	 * @param s Reference to the Simulation object
	 */
	public void applyInitialCondition(Simulation s) {
		totalCells = s.grid.getTotalNumberOfCells();
		numberOfColors = s.getNumberOfColors();
		numberOfComponents = (numberOfColors > 1) ? numberOfColors * numberOfColors - 1 : 1;

		if(numberOfColors != 2) {
			throw new RuntimeException("SU2GlasmaFluxTubes: Invalid number of colors. Only N_c = 2 is supported.");
		}

		transNumCells = GridFunctions.reduceGridPos(s.grid.getNumCells(), direction);
		totalTransCells = GridFunctions.getTotalNumberOfCells(transNumCells);
		aT = s.grid.getLatticeSpacing((direction + 1) % 3);
		g = s.getCouplingConstant();
		this.s = s;

		// Generate transverse gauge fields from nuclei A and B
		SU2GroupElement[][] UA = generateTransverseGaugeFields();
		SU2GroupElement[][] UB = generateTransverseGaugeFields();

		// Compute transverse gauge field and longitudinal electric field from UA, UB
		SU2GroupElement[][] U = new SU2GroupElement[totalTransCells][transNumCells.length];
		for (int i = 0; i < totalTransCells; i++) {
			for (int j = 0; j < transNumCells.length; j++) {
				SU2GroupElement sum = (SU2GroupElement) UA[i][j].add(UB[i][j]);
				U[i][j] = (SU2GroupElement) sum.mult((((SU2GroupElement) sum.adj()).inv()));
			}
		}

		SU2AlgebraElement[] EL =  new SU2AlgebraElement[totalTransCells];
		for (int i = 0; i < totalTransCells; i++) {
			EL[i] = new SU2AlgebraElement();
		}
		SU2GroupElement id = new SU2GroupElement();
		for (int i = 0; i < totalTransCells; i++) {
			SU2GroupElement tmp = new SU2GroupElement();
			for (int j = 0; j < transNumCells.length; j++) {
				int is = GridFunctions.shift(i, j, -1, transNumCells);
				SU2GroupElement Um1 = (SU2GroupElement) U[i][j].adj().sub(id);
				SU2GroupElement diff1 = (SU2GroupElement) UA[i][j].sub(UB[i][j]);
				SU2GroupElement Um2 = (SU2GroupElement) U[is][j].adj().sub(id);
				SU2GroupElement diff2 = (SU2GroupElement) UA[is][j].sub(UB[is][j]);
				tmp.addAssign(diff1.mult(Um1).sub(Um2.mult(diff2)));
			}

			EL[i] = (SU2AlgebraElement) tmp.proj().mult(- 1.0 / (2 * g * aT * aT) * s.grid.getLatticeUnitFactor(direction));
		}

		// Place the fields on the grid
		for (int i = 0; i < totalCells; i++) {
			int[] gridPos = s.grid.getCellPos(i);
			int[] transGridPos = GridFunctions.reduceGridPos(gridPos, direction);
			int transIndex = GridFunctions.getCellIndex(transGridPos, transNumCells);

			int ts = 0;
			for (int j = 0; j < s.getNumberOfDimensions(); j++) {
				if(j != direction) {
					s.grid.setU(i, j, U[transIndex][ts]);
					s.grid.setUnext(i, j, U[transIndex][ts]);
					ts++;
				} else {
					// Manually evolve longitudinal gauge field one time-step.
					s.grid.setUnext(i, direction, EL[transIndex].mult(-s.getTimeStep()).getLink());
				}
			}

			s.grid.setE(i, direction, EL[transIndex]);
		}
	}

	/**
	 * Generates the transverse gauge fields in the temporal gauge based on the MV model.
	 * @return
	 */
	private SU2GroupElement[][] generateTransverseGaugeFields() {
		SU2GroupElement[][] U = new SU2GroupElement[totalTransCells][transNumCells.length];
		SU2GroupElement[] V = generateWilsonLine();

		for (int i = 0; i < totalTransCells; i++) {
			GroupElement V1 = V[i];
			for (int j = 0; j < transNumCells.length; j++) {
				int shifted = GridFunctions.shift(i, j, 1, transNumCells);
				U[i][j] = (SU2GroupElement) V1.mult(V[shifted].adj());
			}
		}

		return U;
	}

	/**
	 * Generates a Wilson line V(x_T) based on the MV model.
	 * @return
	 */
	private SU2GroupElement[] generateWilsonLine() {
		SU2AlgebraElement[] phiAlgebra = new SU2AlgebraElement[totalTransCells];
		for (int i = 0; i < totalTransCells; i++) {
			phiAlgebra[i] = new SU2AlgebraElement();
		}

		for (int c = 0; c < numberOfComponents; c++) {
			double[] rho = generateRandomChargeDensity(mu * g / aT);
			rho = FourierFunctions.regulateChargeDensity2D(rho, transNumCells, UVT, IR, aT);
			double[] phi = FourierFunctions.solvePoisson2D(rho, transNumCells, aT);

			for (int i = 0; i < totalTransCells; i++) {
				phiAlgebra[i].set(c, - phi[i] * g);
			}
		}

		SU2GroupElement[] V = new SU2GroupElement[totalTransCells];
		for (int i = 0; i < totalTransCells; i++) {
			V[i] = (SU2GroupElement) phiAlgebra[i].getLink();
		}

		return V;
	}

	/**
	 * Generates an unregulated charge density distribution for one color component based on the MV-model.
	 * @param gaussianWidth Width of the Gaussian distribution of charges
	 * @return
	 */
	private double[] generateRandomChargeDensity(double gaussianWidth) {
		double[] rho = new double[totalTransCells];
		for (int i = 0; i < totalTransCells; i++) {
			rho[i] = rand.nextGaussian() * gaussianWidth;
		}

		return rho;
	}
}
