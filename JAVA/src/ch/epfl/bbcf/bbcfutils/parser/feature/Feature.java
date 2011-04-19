package ch.epfl.bbcf.bbcfutils.parser.feature;
public abstract class Feature {


	protected String chromosome;
	public abstract String detail();

	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}

	public String getChromosome() {
		return chromosome;
	}

}
