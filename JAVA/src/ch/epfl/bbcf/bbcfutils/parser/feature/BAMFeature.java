package ch.epfl.bbcf.bbcfutils.parser.feature;


//FIXME don't work as expected
public class BAMFeature extends QualitativeFeature{

	public BAMFeature(String chromosome, int start, int end, float score,
			int strand, String name, String attributes) {
		super(chromosome, start, end, score, strand, name, attributes);
		// TODO Auto-generated constructor stub
	}

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


}
