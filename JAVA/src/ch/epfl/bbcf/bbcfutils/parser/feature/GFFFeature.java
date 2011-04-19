package ch.epfl.bbcf.bbcfutils.parser.feature;




public class GFFFeature extends ExtendedQualitativeFeature{


	public GFFFeature(String chromosome, Integer start, Integer end, String name, String id,
			Float score, Integer strand, String type, String attributes) {
		super(chromosome, start, end, score, strand, name,
				attributes, type, id);
	}

	

	
}
