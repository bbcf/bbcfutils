package ch.epfl.bbcf.bbcfutils.parser.feature;

public class QuantitativeFeature extends Feature{


	private Integer start,end;
	private Float score;

	@Override
	public String detail() {
		return this.getClass().getName()+"\n" +
		"chromosome : "+this.chromosome+"\n" +
		"start : "+this.start+"\n" +
		"end : "+this.end+"\n" +
		"score : "+score;
	}

	public void setStart(Integer start) {
		this.start = start;
	}

	public int getStart() {
		return start;
	}

	public void setEnd(Integer end) {
		this.end = end;
	}

	public Integer getEnd() {
		return end;
	}

	public void setScore(Float score) {
		this.score = score;
	}

	public Float getScore() {
		return score;
	}

}
