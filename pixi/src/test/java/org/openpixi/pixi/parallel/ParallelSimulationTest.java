package org.openpixi.pixi.parallel;

import junit.framework.TestCase;
import org.openpixi.pixi.physics.Settings;
import org.openpixi.pixi.physics.Simulation;
import org.openpixi.pixi.physics.solver.Boris;
import org.openpixi.pixi.physics.util.ClassCopier;
import org.openpixi.pixi.util.ResultsComparator;

/**
 * Tests the parallel (threaded) version of openpixi.
 * Multi-threaded and also single-threaded simulations and then compares the results.
 */
public class ParallelSimulationTest extends TestCase {

	public void testParallelSimulation() {
		Settings defaultSettings = new Settings();
		defaultSettings.setNumOfParticles(1000);
		defaultSettings.setIterations(100);
		defaultSettings.setParticleSolver(new Boris());

		Simulation singleThreadedSimulation = new Simulation(defaultSettings);

		Settings multiThreadedSettings = ClassCopier.copy(defaultSettings);
		multiThreadedSettings.setNumOfThreads(10);
		Simulation multiThreadedSimulation = new Simulation(multiThreadedSettings);

		ResultsComparator comparator = new ResultsComparator();
		comparator.compare(
				singleThreadedSimulation.particles, multiThreadedSimulation.particles,
				singleThreadedSimulation.grid, multiThreadedSimulation.grid);
	}
}
