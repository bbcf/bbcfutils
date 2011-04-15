package ch.epfl.bbcf.bbcfutils.parser.feature;
public class BAMFeature extends QualitativeFeature{

	private String readName;
	private int start,stop;

	
	@Override
	public String detail() {
		return "BAMFeature : refName : "+chromosome+
		" readName "+readName+" start "+start+" stop "+stop;
	}

	public void setReadName(String readName) {
		this.readName = readName;
	}

	public String getReadName() {
		return readName;
	}

	public void setRefName(String refName) {
		this.chromosome = refName;
	}

	public String getRefName() {
		return chromosome;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public int getStart() {
		return start;
	}

	public void setStop(int stop) {
		this.stop = stop;
	}

	public int getStop() {
		return stop;
	}


}
