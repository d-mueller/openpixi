package org.openpixi.pixi.physics.initial.CGC;

import org.openpixi.pixi.physics.gauge.DoubleFFTWrapper;
import org.openpixi.pixi.physics.util.GridFunctions;

/**
 * A class with some utility functions for solving the Poisson equation and regulation in momentum space.
 */
public class FourierFunctions {

	/**
	 * Applies hard UV and IR regulation to (one color component of) a 3D charge density.
	 *
	 * @param rho		Original charge density
	 * @param numCells	Grid size
	 * @param UVT		Transverse UV cutoff
	 * @param UVL		Longitudinal UV cutoff
	 * @param IR		IR regulator
	 * @param direction	Longitudinal direction
	 * @return			Regulated charge density
	 */
	public static double[] regulateChargeDensityHard(double[] rho, int[] numCells, double UVT, double UVL, double IR,
													 int direction, double as) {
		int totalCells = numCells[0] * numCells[1] * numCells[2];
		int longitudinalNumCells = numCells[direction];
		int[] transverseNumCells = GridFunctions.reduceGridPos(numCells, direction);

		DoubleFFTWrapper fft = new DoubleFFTWrapper(numCells);
		double[] rhoReg = new double[fft.getFFTArraySize()]; // hardcoded for 3D, I know.
		for (int i = 0; i < totalCells; i++) {
			rhoReg[fft.getFFTArrayIndex(i)] = rho[i];
		}
		fft.complexForward(rhoReg);

		// Remove global charge
		rhoReg[0] = 0.0;
		rhoReg[1] = 0.0;

		// Apply momentum regulation
		for (int i = 0; i < totalCells; i++) {
			int[] gridPos = GridFunctions.getCellPos(i, numCells);
			int longPos = gridPos[direction];
			int[] transPos = GridFunctions.reduceGridPos(gridPos, direction);
			int transIndex = GridFunctions.getCellIndex(transPos, transverseNumCells);

			double kTeff2 = computeEffectiveTransverseMomentumSquared(transIndex, transverseNumCells, as);
			double kT2 = computeTransverseMomentumSquared(transIndex, transverseNumCells, as);
			double kL = Math.abs(computeLatticeMomentum1D(longPos, longitudinalNumCells, as));

			// Apply hard UV regulation
			if(kT2 <=  UVT * UVT && kL <= UVL && kT2 > 0) {
				// Apply 'soft' IR regulation
				double regulator = kTeff2 / (kTeff2 + IR * IR);
				rhoReg[fft.getFFTArrayIndex(i)] *= regulator;
				rhoReg[fft.getFFTArrayIndex(i)+1] *= regulator;
			} else {
				rhoReg[fft.getFFTArrayIndex(i)] = 0.0;
				rhoReg[fft.getFFTArrayIndex(i)+1] = 0.0;
			}
		}

		fft.complexInverse(rhoReg, true);
		for (int i = 0; i < totalCells; i++) {
			rho[i] = rhoReg[fft.getFFTArrayIndex(i)];
		}
		return rho;
	}

	/**
	 * Applies a Gaussian UV regulation in the longitudinal direction as well as hard UV and IR regulation in the
	 * transverse plane to (one color component of) a 3D charge density.
	 *
	 * @param rho			Original charge density
	 * @param numCells		Grid size
	 * @param UVT			Transverse UV cutoff
	 * @param longWidth		Longitudinal Gaussian width in momentum space
	 * @param IR			IR regulator
	 * @param direction		Longitudinal direction
	 * @return				Regulated charge density
	 */
	public static double[] regulateChargeDensityGaussian(double[] rho, int[] numCells, double UVT, double longWidth, double IR,
													 int direction, double as) {
		int totalCells = numCells[0] * numCells[1] * numCells[2];
		int longitudinalNumCells = numCells[direction];
		int[] transverseNumCells = GridFunctions.reduceGridPos(numCells, direction);

		DoubleFFTWrapper fft = new DoubleFFTWrapper(numCells);
		double[] rhoReg = new double[fft.getFFTArraySize()]; // hardcoded for 3D, I know.
		for (int i = 0; i < totalCells; i++) {
			rhoReg[fft.getFFTArrayIndex(i)] = rho[i];
		}
		fft.complexForward(rhoReg);

		// Remove global charge
		rhoReg[0] = 0.0;
		rhoReg[1] = 0.0;

		// Apply momentum regulation
		for (int i = 0; i < totalCells; i++) {
			int[] gridPos = GridFunctions.getCellPos(i, numCells);
			int longPos = gridPos[direction];
			int[] transPos = GridFunctions.reduceGridPos(gridPos, direction);
			int transIndex = GridFunctions.getCellIndex(transPos, transverseNumCells);

			double kTeff2 = computeEffectiveTransverseMomentumSquared(transIndex, transverseNumCells, as);
			double kT2 = computeTransverseMomentumSquared(transIndex, transverseNumCells, as);
			double kL = Math.abs(computeLatticeMomentum1D(longPos, longitudinalNumCells, as));

			// Apply hard UV regulation
			if(kT2 <=  UVT * UVT && kT2 > 0) {
				// Apply 'soft' IR regulation
				double regulator = kTeff2 / (kTeff2 + IR * IR);
				double longRegulator = Math.exp(-0.25*kL*kL*longWidth*longWidth);
				rhoReg[fft.getFFTArrayIndex(i)] *= regulator*longRegulator;
				rhoReg[fft.getFFTArrayIndex(i)+1] *= regulator*longRegulator;
			} else {
				rhoReg[fft.getFFTArrayIndex(i)] = 0.0;
				rhoReg[fft.getFFTArrayIndex(i)+1] = 0.0;
			}
		}

		fft.complexInverse(rhoReg, true);
		for (int i = 0; i < totalCells; i++) {
			rho[i] = rhoReg[fft.getFFTArrayIndex(i)];
		}
		return rho;
	}

	/**
	 * Solves the 2D transverse Poisson equation on the lattice.
	 *
	 * @param rho			2D charge density
	 * @param transNumCells	grid size
	 * @param as			transverse lattice spacing
	 * @return
	 */
	public static double[] solvePoisson2D(double[] rho, int[] transNumCells, double as) {
		int totalCells = transNumCells[0] * transNumCells[1];

		DoubleFFTWrapper fft = new DoubleFFTWrapper(transNumCells);
		double[] phiFFT = new double[fft.getFFTArraySize()];
		for (int i = 0; i < totalCells; i++) {
			phiFFT[fft.getFFTArrayIndex(i)] = rho[i];
		}
		fft.complexForward(phiFFT);

		phiFFT[0] = 0.0;
		phiFFT[1] = 0.0;
		for (int i = 1; i < totalCells; i++) {
			double kTeff2 = computeEffectiveTransverseMomentumSquared(i, transNumCells, as);
			phiFFT[fft.getFFTArrayIndex(i)]   *= 1.0 / kTeff2;
			phiFFT[fft.getFFTArrayIndex(i)+1] *= 1.0 / kTeff2;
		}
		fft.complexInverse(phiFFT, true);

		double[] phi = new double[totalCells];
		for (int i = 0; i < totalCells; i++) {
			phi[i] = phiFFT[fft.getFFTArrayIndex(i)];
		}

		return phi;
	}

	/**
	 * Computes the square of the effective transverse momentum on the lattice as used in the Poisson equation.
	 *
	 * @param cellIndex		transverse lattice index
	 * @param transNumCells	transverse grid size
	 * @param as			transverse lattice spacing
	 * @return				effective kT squared
	 */
	public static double computeEffectiveTransverseMomentumSquared(int cellIndex, int[] transNumCells, double as) {
		int[] transversalGridPos = GridFunctions.getCellPos(cellIndex, transNumCells);

		double momentumSquared = 2.0 * transNumCells.length;
		for (int i = 0; i < transNumCells.length; i++) {
			momentumSquared -= 2.0 * Math.cos((2.0 * Math.PI * transversalGridPos[i]) / transNumCells[i]);
		}

		return momentumSquared / (as * as);
	}

	/**
	 * Compute the square of the 'real' transverse momentum on the lattice.
	 *
	 * @param cellIndex		transverse lattice index
	 * @param transNumCells	transverse grid size
	 * @param as			transverse lattice spacing
	 * @return				'real' kT squared
	 */
	public static double computeTransverseMomentumSquared(int cellIndex, int[] transNumCells, double as) {
		int[] transversalGridPos = GridFunctions.getCellPos(cellIndex, transNumCells);
		double twopi = 2.0 * Math.PI;

		double momentumSquared = 0.0;
		for (int i = 0; i < transNumCells.length; i++) {
			double momentumComponent;
			int n = transNumCells[i];
			if(transversalGridPos[i] < n / 2) {
				momentumComponent = twopi * transversalGridPos[i] / (as * n);
			} else {
				momentumComponent = twopi * (n - transversalGridPos[i]) / (as * n);
			}

			momentumSquared += momentumComponent * momentumComponent;
		}

		return momentumSquared;
	}

	/**
	 * Computes the longitudinal momentum on the lattice.
	 *
	 * @param pos		longitudinal lattice position
	 * @param numCells	number of cells in the longitudinal direction
	 * @param as		longitudinal lattice spacing
	 * @return			kL
	 */
	public static double computeLatticeMomentum1D(int pos, int numCells, double as) {
		double delta = pos / ((double) numCells);
		if(delta < 0.5)
		{
			return  2.0 * delta * Math.PI / as;
		}
		return 2.0 * (delta - 1.0) * Math.PI / as;
	}
}
