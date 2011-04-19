package ch.epfl.bbcf.bbcfutils.parser.feature;




public class JSONFeature {

	private String chromosome,id,name,type;
	private Integer start,end,strand;
	private Float score;



	public JSONFeature(QualitativeFeature feature) {
		setChromosome(feature.getChromosome());
		setName(feature.getName());
		setStart(feature.getStart());
		setEnd(feature.getEnd());
		setStrand(feature.getStrand());
		setScore(feature.getScore());
	}
	
	public JSONFeature(ExtendedQualitativeFeature feature) {
		setChromosome(feature.getChromosome());
		setName(feature.getName());
		setStart(feature.getStart());
		setEnd(feature.getEnd());
		setStrand(feature.getStrand());
		setType(feature.getType());
		setId(feature.getIdentifier());
		setScore(feature.getScore());
	}
	public void setId(String id) {
		this.id = id;
	}
	public String getId() {
		return id;
	}
	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}
	public String getChromosome() {
		return chromosome;
	}
	public void setName(String name) {
		this.name = name;
	}
	public String getName() {
		return name;
	}
	public void setType(String type) {
		this.type = type;
	}
	public String getType() {
		return type;
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
	public void setStrand(Integer strand) {
		this.strand = strand;
	}
	public Integer getStrand() {
		return strand;
	}

	public void setScore(Float score) {
		this.score = score;
	}

	public Float getScore() {
		return score;
	}



}
