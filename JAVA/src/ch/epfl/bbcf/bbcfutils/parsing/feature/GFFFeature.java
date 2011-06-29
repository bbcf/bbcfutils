package ch.epfl.bbcf.bbcfutils.parsing.feature;




public class GFFFeature extends Feature{


	private String chromosome;
	private Integer start;
	private Integer end;
	private String name;
	private String id;
	private Float score;
	private Integer strand;
	private String type;
	private String attributes;
	private String source;
	private Integer frame;

	

	public GFFFeature(String chromosome, String source,Integer start, Integer end, String name, String id,
			Float score, Integer strand, String type, String attributes) {
		this.chromosome= chromosome;
		this.source = source;
		this.start = start;
		this.end = end;
		this.name = name;
		this.id = id;
		this.score = score;
		this.strand = strand;
		this.type = type;
		this.attributes = attributes;
	}

	@Override
	public
	String display() {
		return ds(this.chromosome)+SEP+ds(this.source)+SEP+ds(this.type)+ds(this.start)
		+SEP+ds(this.end)+SEP+ds(this.score)+SEP+dss(this.strand)+SEP+ds(this.frame)
		+ds(attributes);
	}

	
	public String getChromosome() {
		return chromosome;
	}

	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
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

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public Float getScore() {
		return score;
	}

	public void setScore(Float score) {
		this.score = score;
	}

	public Integer getStrand() {
		return strand;
	}

	public void setStrand(Integer strand) {
		this.strand = strand;
	}

	public String getType() {
		return type;
	}

	public void setType(String type) {
		this.type = type;
	}

	public String getAttributes() {
		return attributes;
	}

	public void setAttributes(String attributes) {
		this.attributes = attributes;
	}
	
	
	
}
