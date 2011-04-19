package ch.epfl.bbcf.bbcfutils.parser.feature;

public class ExtendedQualitativeFeature extends QualitativeFeature{

	private String type,identifier;
	public ExtendedQualitativeFeature(String chromosome, Integer start, Integer end,
			Float score, Integer strand, String name, String attributes,String type,String id) {
		super(chromosome, start, end, score, strand, name, attributes);
		this.type=type;
		this.setIdentifier(id);
		
	}
	public ExtendedQualitativeFeature() {
	}
	public void setType(String type) {
		this.type = type;
	}
	public String getType() {
		return type;
	}
	public void setIdentifier(String identifier) {
		this.identifier = identifier;
	}
	public String getIdentifier() {
		return identifier;
	}

}
