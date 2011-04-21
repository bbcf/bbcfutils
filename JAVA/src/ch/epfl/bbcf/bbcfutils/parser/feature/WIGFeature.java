package ch.epfl.bbcf.bbcfutils.parser.feature;

public class WIGFeature extends QuantitativeFeature{

	public static final int VARIABLE_STEP = 1;
	public static final int FIXED_STEP = 2;
	private int format;
	
	
	private Integer step,span;
	
	
	
	public WIGFeature(WIGFeature feat) {
		super(feat.chromosome,feat.getStart(),feat.getEnd(),feat.getScore());
		this.step=feat.step;
		this.span=feat.span;
		this.format=feat.format;
	}


	public WIGFeature(){}

	public WIGFeature(String chr, Integer start2, Integer end2, Float score2) {
		super(chr,start2,end2,score2);
	}




	public void setFormat(int format) {
		this.format = format;
	}
	
	@Override
	public String detail(){
		return super.detail()+" step : "+this.step+" format : "+this.format+
		" span : "+this.span;
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
	public int getSpan() {
		return span;
	}
	public void setStep(int step) {
		this.step = step;
	}
	public int getStep() {
		return step;
	}
	

	
}
