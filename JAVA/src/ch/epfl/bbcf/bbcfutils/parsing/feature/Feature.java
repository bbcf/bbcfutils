package ch.epfl.bbcf.bbcfutils.parsing.feature;

public abstract class Feature {

	protected static final String SEP = "\t";
	public abstract String display();
	public abstract String getChromosome();
	
	/**
	 * return the display of the parameter
	 * @param str
	 * @return
	 */
	protected String ds(String str){
		if(null==str){
			return ".";
		}
		return str;
	}
	
	

	/**
	 * return the display of the parameter
	 * @param str
	 * @return
	 */
	protected String ds(Integer str){
		if(null==str){
			return ".";
		}
		return str.toString();
	}
	

	/**
	 * return the display of the parameter
	 * @param str
	 * @return
	 */
	protected String ds(Float str){
		if(null==str){
			return ".";
		}
		return str.toString();
	}
	
	/**
	 * display the strand
	 * @param strand
	 * @return
	 */
	protected String dss(Integer strand){
		if(strand==null || strand==0){
			return ".";
		}
		if(strand>0){
			return "+";
		} else {
			return "-";
		}
	}

}
