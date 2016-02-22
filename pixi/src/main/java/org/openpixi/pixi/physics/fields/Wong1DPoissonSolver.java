package org.openpixi.pixi.physics.fields;

import org.openpixi.pixi.math.ElementFactory;
import org.openpixi.pixi.physics.gauge.DoubleFFTWrapper;
import org.openpixi.pixi.physics.grid.Grid;

public class Wong1DPoissonSolver implements PoissonSolver {

	int numberOfCells;
	double as;

	public void solve (Grid g) {
		numberOfCells = g.getTotalNumberOfCells();
		as = g.getLatticeSpacing();
		ElementFactory factory = g.getElementFactory();
		DoubleFFTWrapper fft = new DoubleFFTWrapper(g.getNumCells());
		// Assume that gauge fields are zero initially. The Gauss constraint then reduces to the abelian case.
		for (int j = 0; j < factory.numberOfComponents; j++) {
			// Fourier transform charge density and solve in momentum space.
			double[] fftArray = new double[fft.getFFTArraySize()];
			for (int i = 0; i < numberOfCells; i++) {
				fftArray[fft.getFFTArrayIndex(i)] = g.getRho(i).get(j);
			}
			fft.complexForward(fftArray);

			for (int i = 1; i < numberOfCells; i++) {
				double invLaplace = 1.0 / computeLatticeMomentumSquared(i);
				fftArray[fft.getFFTArrayIndex(i)] *= invLaplace;
				fftArray[fft.getFFTArrayIndex(i) + 1] *= invLaplace;
			}
			fftArray[0] = 0.0;
			fftArray[1] = 0.0;

			// Inverse Fourier transform.
			fft.complexInverse(fftArray, true);

			// Compute chromoelectric field from potential.
			for (int i = 1; i < numberOfCells; i++) {
				int l = g.shift(i, 0, -1);
				int r = g.shift(i, 0, 1);
				double E = (fftArray[fft.getFFTArrayIndex(r)] - fftArray[fft.getFFTArrayIndex(l)]) / (2 * as);
				g.getE(i, 0).set(j, E);
			}
		}

	}

	private double computeLatticeMomentumSquared(int cellIndex) {
		return 2.0 * (1.0 - Math.cos((2.0 * Math.PI * cellIndex) / numberOfCells)) / (as * as);
	}

}
