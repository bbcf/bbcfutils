package ch.epfl.bbcf.bbcfutils.parsing.feature;

import java.util.HashMap;
import java.util.Map;

public class Track {

	private Map<String,String> attributes;



	public Track(){
		this.attributes = new HashMap<String, String>();
	}
	
	protected void setAttributes(Map<String,String> attributes) {
		this.attributes = attributes;
	}
	
	public void addAttribute(String key,String value) {
		this.attributes.put(key,value);
	}
	/**
	 * get the attributes of the track
	 * @return an HASHMap containing attributes of the line
	 */
	public Map<String,String> getAttributes() {
		return attributes;
	}

	public String detail() {
		String detail="";
		for(Map.Entry<String, String>entry : attributes.entrySet()){
			detail+=entry.getKey()+"="+entry.getValue()+" ";
		}
		return detail;
	}
}
