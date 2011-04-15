package ch.epfl.bbcf.bbcfutils.parser.feature;

import java.util.HashMap;
import java.util.Map;

public class QualitativeFeature extends Feature{

	protected String attributes,name;
	protected int start,end,strand;
	protected float score;

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

	public void setScore(float score) {
		this.score = score;
	}

	public float getScore() {
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
