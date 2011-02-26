package ch.epfl.bbcf.feature;


public class WIGFeature implements Feature{

	public static final int VARIABLE_STEP = 1;
	public static final int FIXED_STEP = 2;
	private int format;
	
	private String chromosome;
	
	private int start,end,step,span;
	
	private float score;
	
	
	public WIGFeature(WIGFeature cur_feature) {
		this.start=cur_feature.start;
		this.end=cur_feature.end;
		this.score=cur_feature.score;
		this.chromosome=cur_feature.chromosome;
		this.step=cur_feature.step;
		this.span=cur_feature.span;
		this.format=cur_feature.format;
	}

	public WIGFeature() {
	}

	public void setFormat(int format) {
		this.format = format;
	}
	
	@Override
	public String detail(){
		return "WIGFeature : chr : "+this.chromosome+
		" start : "+this.start+" end : "+this.end+
		" step : "+this.step+" score : "+this.score+
		" span : "+this.span;
	}
	
	/**
	 * get the format of the feature : variable step or fixed step
	 * @return 1 for variable, 2 for fixed
	 */
	public int getFormat() {
		return format;
	}
	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}
	public String getChromosome() {
		return chromosome;
	}
	public void setStart(int start) {
		this.start = start;
	}
	public int getStart() {
		return start;
	}
	public void setSpan(int span) {
		this.span = span;
	}
	public int getSpan() {
		return span;
	}
	public void setStep(int step) {
		this.step = step;
	}
	public int getStep() {
		return step;
	}
	public void setEnd(int end) {
		this.end = end;
	}
	public int getEnd() {
		return end;
	}
	public void setScore(float score) {
		this.score = score;
	}
	public float getScore() {
		return score;
	}

	
}
