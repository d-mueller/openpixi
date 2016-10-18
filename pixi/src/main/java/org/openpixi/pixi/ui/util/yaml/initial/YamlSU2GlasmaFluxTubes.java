package org.openpixi.pixi.ui.util.yaml.initial;

import org.openpixi.pixi.physics.initial.IInitialCondition;
import org.openpixi.pixi.physics.initial.YM.SU2GlasmaFluxTubes;

public class YamlSU2GlasmaFluxTubes {

	/**
	 * Longitudinal direction of the glasma fields
	 */
	public Integer direction;

	/**
	 * MV-Parameter
	 */
	public Double mu;

	/**
	 * Transverse UV cutoff
	 */
	public Double ultravioletCutoffTransverse;

	/**
	 * Infrared regulator in the Poisson equation
	 */
	public Double infraredCoefficient;

	/**
	 * Optional seed for the random number generator
	 */
	public Integer randomSeed;

	public IInitialCondition getInitialCondtion() {
		boolean useSeed = (randomSeed != null);
		if(!useSeed) {
			randomSeed = 0;
		}

		return new SU2GlasmaFluxTubes(direction, mu, infraredCoefficient, ultravioletCutoffTransverse, useSeed, randomSeed);
	}
}
