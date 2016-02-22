package org.openpixi.pixi.ui.util.yaml;

import java.awt.Color;
import java.lang.reflect.Field;
import java.util.List;

import org.openpixi.pixi.physics.Settings;
import org.openpixi.pixi.physics.particles.Wong1DParticle;
import org.openpixi.pixi.physics.particles.YangMillsParticle;

public class YamlYangMillsParticle {
	public List<Double> position;
	public List<Double> velocity;
	public Double r;
	public Double m;
	public List<Double> q;
	public String color;

	public void applyTo(Settings settings) {

		Wong1DParticle p = getParticle(settings);

		settings.addParticle(p);
	}

	/**
	 * Creates a new particle and applies the settings from the
	 * YAML document to it.
	 * @return new particle
	 */
	public Wong1DParticle getParticle(Settings settings) {
		Wong1DParticle p = new Wong1DParticle(settings.getNumberOfDimensions(), settings.getNumberOfColors());


		if (position != null) {
			for (int i = 0; i < position.size(); i++) {
				p.setPosition(i, position.get(i));
				p.setPrevPosition(i, position.get(i));
			}
		}

		if(velocity != null)
			for(int i = 0; i < velocity.size(); i++)
				p.setVelocity(i, velocity.get(i));

		if (r != null) {
			p.setRadius(r);
		}

		if (m != null) {
			p.mass = m;
		}

		double unitFactor = settings.getCouplingConstant() * settings.getGridStep();
		if (q != null)
			for(int c = 0; c < q.size(); c++) {
				p.Q0.set(c, q.get(c) * unitFactor);
				p.Q1.set(c, q.get(c) * unitFactor);
			}

		if (color != null) {
			p.setDisplayColor(getColorFromString(color));
		}
		return p;
	}

	/**
	 * Convert a color string from the YAML file into Java color.
	 * @param colorstring This can be a name (e.g. "red", "blue", ...),
	 * a HEX code (e.g. "FFFFFF"), or "random".
	 * @return
	 */
	static private Color getColorFromString(String colorstring) {
		// Default is black
		Color c = null;

		// Check if it is a Java string (e.g. "blue", "red", ...)
		try {
			// Access Color.blue, Color.red, ... by reflection:
			Field field = Color.class.getField(colorstring);
			c = (Color)field.get(null);
		} catch (NoSuchFieldException e) {
		} catch (IllegalAccessException e) {
		}

		if (c == null) {
			// Check if it is a HEX color code (e.g. "FFFFFF")
			try {
				c = Color.decode("#" + colorstring);
			} catch (NumberFormatException e) {
			}
		}

		if (colorstring.equals("random")) {
			// Use random color
			c = new Color((int)(Math.random() * 0x1000000));
		}

		if (c == null) {
			// Fallback color
			c = Color.black;
		}
		return c;
	}

}
