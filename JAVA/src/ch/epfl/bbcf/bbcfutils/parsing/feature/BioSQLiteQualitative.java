package ch.epfl.bbcf.bbcfutils.parsing.feature;

import java.util.HashMap;
import java.util.Map;

public class BioSQLiteQualitative extends Feature{

	protected String attributes,name,chromosome;
	protected Integer start,end,strand;
	protected Float score;


	public BioSQLiteQualitative(String chromosome,
			Integer start,Integer end,Float score,Integer strand,
			String name,String attributes){
		this.chromosome = chromosome;
		this.start = start;
		this.end = end;
		this.score = score;
		this.strand = strand;
		this.name = name;
		this.attributes = attributes;
	}

	public BioSQLiteQualitative() {
	}

	public static String buildAttributesfromBEDFeature(
			Integer thickStart, Integer thickEnd, 
			String itemRgb, String blockCount,
			String blockSizes, String blockStarts){
		String attributes ="";
		if(thickStart!=null){
			attributes+=formatAttribute("thickStart", thickStart);
		}
		if(thickEnd!=null){
			attributes+=formatAttribute("thickEnd", thickEnd);
		}
		if(itemRgb!=null){
			attributes+=formatAttribute("itemRgb", itemRgb);
		}
		if(blockCount!=null){
			attributes+=formatAttribute("blockCount",blockCount);
		}
		if(blockSizes!=null){
			attributes+=formatAttribute("blockSizes",blockSizes);
		}
		if(blockStarts!=null){
			attributes+=formatAttribute("blockStarts",blockStarts);
		}
		return attributes;
	}



	public BEDFeature toBEDFeature() {
		Map<String,String> atts = getMapAttributes();
		Integer ths = null;
		Integer the = null;
		
		
		
		try {
			ths = Integer.parseInt(atts.get("thickStart"));
		} catch (NumberFormatException e){};
		try {
			the = Integer.parseInt(atts.get("thickEnd"));
		} catch (NumberFormatException e){};

		BEDFeature bed = new BEDFeature(this.chromosome, this.start, this.end, this.name, this.strand, this.score, 
				ths, the,
				atts.get("itemRgb"),atts.get("blockCount"), 
				atts.get("blockSizes"),atts.get("blockStarts"));
		return bed;
	}




	private static String formatAttribute(String name,int attribute){
		return name+" \""+attribute+"\";";
	}
	private static String formatAttribute(String name,String attribute){
		return name+" \""+attribute+"\";";
	}

	public void setAttributes(String attributes) {
		this.attributes = attributes;
	}

	public String getAttributes() {
		return attributes;
	}

	public Map<String,String> getMapAttributes(){
		Map<String, String> attMap = new HashMap<String, String>();
		String[]atts = attributes.split(";");
		for(String att:atts){
			String[] keyvalue = att.split("\\s");
			if(keyvalue.length>1){
				String str = keyvalue[1];
				str = str.replaceAll("\"", "");
				attMap.put(keyvalue[0],str);
			}
		}
		return attMap;

	}

	public void setName(String name) {
		this.name = name;
	}

	public String getName() {
		return name;
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
	public String getChromosome(){
		return this.chromosome;
	}
	public void setChromosome(String chr){
		this.chromosome=chr;
	}
	@Override
	public String display() {
		return this.chromosome+SEP+this.start+SEP+this.end+
		SEP+this.name+SEP+this.score+
		SEP+this.strand+
		SEP+this.attributes;
	}




}
