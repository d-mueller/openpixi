package org.openpixi.pixi.physics.force;

import org.openpixi.pixi.physics.Simulation;
import org.openpixi.pixi.physics.particles.IParticle;
import org.openpixi.pixi.physics.particles.Wong1DParticle;


public class Wong1DForce implements Force {


    private int numberOfDimensions;
    private int numberOfColors;
    private int numberOfComponents;

    private double couplingConstant;

    public Wong1DForce(Simulation s)
    {
        this(s.getNumberOfDimensions(), s.getNumberOfColors(), s.getCouplingConstant());
    }

    public Wong1DForce(int numberOfDimensions, int numberOfColors, double couplingConstant)
    {
        this.numberOfDimensions = numberOfDimensions;
        this.numberOfColors = numberOfColors;
        this.couplingConstant = couplingConstant;

        if(this.numberOfColors > 1)
        {
            this.numberOfComponents = this.numberOfColors * this.numberOfColors - 1;
        }
        else
        {
            this.numberOfComponents = 1;
        }
    }

    public Wong1DForce()
    {
        this(1, 1, 1.0);
    }

    public double getForce(int i, IParticle p)
    {
        Wong1DParticle P = (Wong1DParticle) p;
        double f = 0.0;
        for(int c = 0; c < this.numberOfComponents; c++)
        {
            f += P.Q0.get(c) * P.E.get(c);
        }
        f *= this.couplingConstant;
        return f;
    }

}
