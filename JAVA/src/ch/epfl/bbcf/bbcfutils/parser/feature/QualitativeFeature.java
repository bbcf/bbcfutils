package ch.epfl.bbcf.bbcfutils.parser.feature;

import java.util.HashMap;
import java.util.Map;

public class QualitativeFeature extends Feature{

	protected String attributes,name;
	protected Integer start,end,strand;
	protected Float score;


	public QualitativeFeature(String chromosome,
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

	public QualitativeFeature() {
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
				attMap.put(keyvalue[0],keyvalue[1]);
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

	@Override
	public String detail() {
		return this.getClass().getName()+" : chr : "+this.chromosome+
		" start : "+this.start+" end : "+this.end+
		" name : "+this.name+" score : "+this.score+
		" strand : "+this.strand+
		" attributes : "+this.attributes;
	}







}
