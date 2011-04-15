package ch.epfl.bbcf.bbcfutils.parser.feature;




public class JSONFeature {

	private String chromosome,id,name,type;
	private int start,end,strand;



	public JSONFeature(QualitativeFeature feature) {
		setChromosome(feature.getChromosome());
		setName(feature.getName());
		setStart(feature.getStart());
		setEnd(feature.getEnd());
		setStrand(feature.getStrand());
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
	public void setStart(int start) {
		this.start = start;
	}
	public int getStart() {
		return start;
	}
	public void setEnd(int end) {
		this.end = end;
	}
	public int getEnd() {
		return end;
	}
	public void setStrand(int strand) {
		this.strand = strand;
	}
	public int getStrand() {
		return strand;
	}



}
