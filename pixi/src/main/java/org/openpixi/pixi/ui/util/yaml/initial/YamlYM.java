package org.openpixi.pixi.ui.util.yaml.initial;

import org.openpixi.pixi.physics.Settings;

import java.util.ArrayList;

public class YamlYM {


	public ArrayList<YamlSU2PlaneWave> SU2PlaneWaves = new ArrayList<YamlSU2PlaneWave>();
	public ArrayList<YamlSU2GlasmaFluxTubes> SU2GlasmaFluxTubes = new ArrayList<YamlSU2GlasmaFluxTubes>();

	public void applyTo(Settings s) {
		for (YamlSU2PlaneWave wave : SU2PlaneWaves) {
			if (wave.checkConsistency(s)) {
				s.addInitialConditions(wave.getFieldGenerator());
			}
		}

		for (YamlSU2GlasmaFluxTubes glasma : SU2GlasmaFluxTubes) {
			s.addInitialConditions(glasma.getInitialCondtion());
		}
	}
}
