package ch.epfl.bbcf.bbcfutils.parsing.feature;

public class BioSQLiteQualitativeExt extends BioSQLiteQualitative{

	private String type,identifier;
	
	
	
	public BioSQLiteQualitativeExt(String chromosome, Integer start, Integer end,
			Float score, Integer strand, String name, String attributes,String type,String id) {
		super(chromosome, start, end, score, strand, name, attributes);
		this.type=type;
		this.setIdentifier(id);

	}
	public BioSQLiteQualitativeExt() {
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
	public String getChromosome(){
		return this.chromosome;
	}

	public GFFFeature toGFFFeature(){
		return new GFFFeature(chromosome,"GDV", start, end, name, identifier, score, strand, type, attributes);
	}

}
