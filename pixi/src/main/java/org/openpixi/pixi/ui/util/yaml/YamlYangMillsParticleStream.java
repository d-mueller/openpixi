package org.openpixi.pixi.ui.util.yaml;

import java.util.List;
import java.util.Random;

import org.openpixi.pixi.physics.Settings;
import org.openpixi.pixi.physics.particles.IParticle;

public class YamlYangMillsParticleStream {
	public YamlYangMillsParticle particle;
	public List<Double> distances;
	public List<Double> randomPositions;
	public List<Double> randomVelocities;
	public List<Double> randomGaussPositions;
	public List<Double> randomGaussVelocities;
	public Integer number;

	private int numberOfDimensions;
	private Random random = new Random();


	/**
	 * Creates a stream of particles. The particle is copied
	 * 'number' times. Each time, the positions are
	 * adjusted by adding the distances defined in List<Double> distances.
	 * @param settings Settings object to which particles are added.
	 */
	public void applyTo(Settings settings) {
		IParticle p;
		numberOfDimensions = settings.getNumberOfDimensions();
		double[] distances = new double[numberOfDimensions];
		double number = 0;

		if(this.distances != null)
			for(int i = 0; i < numberOfDimensions; i++)
				distances[i] = this.distances.get(i);

		if (this.number != null) {
			number = this.number;
		}

		if (this.particle != null) {

			for (int n = 0; n < number; n++)
			{
				p = particle.getParticle(settings);
				for(int i = 0; i < numberOfDimensions; i++) {
					p.addPosition(i, n * distances[i]);
					p.addPrevPosition(i, n * distances[i]);
				}

				applyRandomModifications(p);
				settings.addParticle(p);
			}
		}
	}

	private void applyRandomModifications(IParticle p) {

		if(this.randomPositions != null)
			for(int i = 0; i < numberOfDimensions; i++) {
				double r = random.nextDouble() * this.randomPositions.get(i);
				p.addPosition(i, r);
				p.addPrevPosition(i, r);
			}

		if(this.randomVelocities != null)
			for(int i = 0; i < numberOfDimensions; i++)
				p.addVelocity(i, random.nextDouble() * this.randomVelocities.get(i));

		if(this.randomGaussPositions != null)
			for(int i = 0; i < numberOfDimensions; i++) {
				double r = random.nextGaussian() * this.randomGaussPositions.get(i);
				p.addPosition(i, r);
				p.addPrevPosition(i, r);
			}

		if(this.randomGaussVelocities != null)
			for(int i = 0; i < numberOfDimensions; i++)
				p.addVelocity(i, random.nextGaussian() * this.randomGaussVelocities.get(i));
	}
}
