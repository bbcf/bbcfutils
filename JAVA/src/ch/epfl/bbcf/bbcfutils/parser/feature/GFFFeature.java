package ch.epfl.bbcf.bbcfutils.parser.feature;




public class GFFFeature extends QualitativeFeature{

	private String id;
	private String type;
	@SuppressWarnings("unused")
	private String attributes;

	public GFFFeature(String chromosome, int start, int end, String name, String id,
			float score, int strand, String type, String attributes) {
		this.setChromosome(chromosome);
		this.setStart(start);
		this.setEnd(end);
		this.setName(name);
		this.setId(id);
		this.setScore(score);
		this.setStrand(strand);
		this.setType(type);
		this.setAttributes(attributes);

	}

	@Override
	public String detail() {
		return super.detail()+"\n"+
		"type "+type+" \n" +
		"id "+id;
	}

	public void setId(String id) {
		this.id = id;
	}

	public String getId() {
		return id;
	}

	public void setType(String type) {
		this.type = type;
	}

	public String getType() {
		return type;
	}
}
