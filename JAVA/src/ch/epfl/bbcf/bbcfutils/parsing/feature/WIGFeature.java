package ch.epfl.bbcf.bbcfutils.parsing.feature;

public class WIGFeature extends Feature{

	public static final int VARIABLE_STEP = 1;
	public static final int FIXED_STEP = 2;
	private int format;
	
	
	private String chromosome;
	

	private Integer step,span,start,end;
	

	

	public WIGFeature(WIGFeature feat) {
		this.chromosome=feat.chromosome;
		this.start=feat.start;
		this.end = feat.end;
		this.score=feat.score;
		this.step=feat.step;
		this.span=feat.span;
		this.format=feat.format;
	}

	public WIGFeature(){}

	public WIGFeature(String chr, Integer start2, Integer end2, Float score2) {
		this.chromosome=chr;
		this.start=start2;
		this.end = end2;
		this.score=score2;
	}




	public void setFormat(int format) {
		this.format = format;
	}
	
	
	/**
	 * get the format of the feature : variable step or fixed step
	 * @return 1 for variable, 2 for fixed
	 */
	public int getFormat() {
		return format;
	}

	public void setSpan(int span) {
		this.span = span;
	}
	public Integer getSpan() {
		return span;
	}
	public void setStep(int step) {
		this.step = step;
	}
	public int getStep() {
		return step;
	}

	@Override
	public String display() {
		return this.chromosome+SEP+this.start+SEP+this.end+SEP+this.score;
	}

	public void setScore(float score2) {
		this.score = score2;
	}
	
	public String getChromosome() {
		return chromosome;
	}

	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}
	public void setScore(Float score) {
		this.score = score;
	}
	
	public Integer getStart() {
		return start;
	}



	public void setStart(Integer start) {
		this.start = start;
	}



	public Integer getEnd() {
		return end;
	}



	public void setEnd(Integer end) {
		this.end = end;
	}
	private Float score;
	
	
	public Float getScore() {
		return score;
	}
}
