package ch.epfl.bbcf.bbcfutils.parsing.feature;


//FIXME don't work as expected
public class BAMFeature extends Feature{

	public BAMFeature(String chromosome, int start, int end, float score,
			int strand, String name, String attributes) {
	}

	private String readName,chromosome;
	private int start,stop;

	


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

	public void setStart(Integer start) {
		this.start = start;
	}

	public Integer getStart() {
		return start;
	}

	public void setStop(Integer stop) {
		this.stop = stop;
	}

	public int getStop() {
		return stop;
	}

	@Override
	public String display() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getChromosome() {
		// TODO Auto-generated method stub
		return null;
	}


}
