package ch.epfl.bbcf.bbcfutils.parsing.feature;

public class BioSQLiteQuantitative extends Feature{

	private String chromosome;
	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}


	private Integer start,end;
	private Float score;

	public BioSQLiteQuantitative(String chromosome,Integer start,Integer end,Float score){
		this.chromosome=chromosome;
		this.start=start;
		this.end=end;
		this.score=score;
	}
	public BioSQLiteQuantitative() {
	}
	

	public void setStart(Integer start) {
		this.start = start;
	}

	public Integer getStart() {
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
	public String getChromosome(){
		return this.chromosome;
	}
	@Override
	public String display() {
		return this.chromosome+SEP+this.start+SEP+this.end+SEP+score;
	}
	
	
	public WIGFeature toWIGFeature(){
		return new WIGFeature(chromosome, start, end, score);
	}

}
