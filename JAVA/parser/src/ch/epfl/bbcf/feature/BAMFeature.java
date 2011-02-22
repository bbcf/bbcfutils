package ch.epfl.bbcf.feature;

public class BAMFeature implements Feature{

	private String refName,readName;
	private int start,stop;

	
	@Override
	public String detail() {
		return "BAMFeature : refName : "+refName+
		" readName "+readName+" start "+start+" stop "+stop;
	}

	public void setReadName(String readName) {
		this.readName = readName;
	}

	public String getReadName() {
		return readName;
	}

	public void setRefName(String refName) {
		this.refName = refName;
	}

	public String getRefName() {
		return refName;
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

	@Override
	public String getChromosome() {
		return refName;
	}

	@Override
	public void setChromosome(String chromosome) {
		this.refName = chromosome;
	}

	

}
