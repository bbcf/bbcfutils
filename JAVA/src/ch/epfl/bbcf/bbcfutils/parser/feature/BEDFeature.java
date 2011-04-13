package ch.epfl.bbcf.bbcfutils.parser.feature;

public class BEDFeature implements Feature{

	private String chromosome;
	private int start,end,thickStart,thickEnd;
	private String name,itemRgb,blockCount,blockSize,blockStarts;
	private Float score;
	private Integer strand;
	
	public BEDFeature(String chromosome2, int start2, int end2, String name2,
			Integer strand2, Float score2, int thickStart, int thickEnd, 
			String itemRgb, String blockCount, String blockSizes, String blockStarts) {
		this.chromosome = chromosome2;
		this.start = start2;
		this.end = end2;
		this.name = name2;
		this.strand = strand2;
		this.score = score2;
		this.thickStart=thickStart;
		this.thickEnd=thickEnd;
		this.itemRgb=itemRgb;
		this.blockCount=blockCount;
		this.blockSize=blockSizes;
		this.blockStarts=blockStarts;
	}
	
	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}
	/**
	 * get the chromosome of the feature
	 * @return the chromosome
	 */
	public String getChromosome() {
		return chromosome;
	}
	/**
	 * set the start
	 * @param start
	 */
	public void setStart(int start) {
		this.start = start;
	}
	/**
	 * get the start base of this feature
	 * @return the start
	 */
	public int getStart() {
		return start;
	}
	public void setName(String name) {
		this.name = name;
	}
	/**
	 * get the name of this feature - can be null
	 * @return the name
	 */
	public String getName() {
		return name;
	}
	/**
	 * get the end point of this feature
	 * @param end
	 */
	public void setEnd(int end) {
		this.end = end;
	}
	public int getEnd() {
		return end;
	}
	
	public void setScore(float score) {
		this.score = score;
	}
	/**
	 * get the score of the feature - can be null
	 * @return the score
	 */
	public Float getScore() {
		if(null==score){
			return 0f;
		}
		return score;
	}
	public void setStrand(int strand) {
		this.strand = strand;
	}
	/**
	 * get the stand 1 or -1 for this feature - can be null
	 * @return the strand
	 */
	public int getStrand() {
		if(null==strand){
			return 0;
		}
		return strand;
	}
	@Override
	public String detail(){
		return "BEDFeature : chr : "+this.chromosome+
		" start : "+this.start+" end : "+this.end+
		" name : "+this.name+" score : "+this.score+
		" strand : "+this.strand+
		" thickStart : "+this.thickStart+
		" thickEnd : "+this.thickEnd+
		" itemRgb : "+this.itemRgb+
		" blockCount : "+this.blockCount+
		" blockSize : "+this.blockSize+
		" blockStarts : "+this.blockStarts;
	}

	public void setThickStart(int thickStart) {
		this.thickStart = thickStart;
	}

	public int getThickStart() {
		return thickStart;
	}

	public void setThickEnd(int thickEnd) {
		this.thickEnd = thickEnd;
	}

	public int getThickEnd() {
		return thickEnd;
	}

	public void setBlockCount(String blockCount) {
		this.blockCount = blockCount;
	}

	public String getBlockCount() {
		return blockCount;
	}

	public void setBlockSize(String blockSize) {
		this.blockSize = blockSize;
	}

	public String getBlockSize() {
		return blockSize;
	}

	public void setItemRgb(String itemRgb) {
		this.itemRgb = itemRgb;
	}

	public String getItemRgb() {
		return itemRgb;
	}

	public void setBlockStarts(String blockStarts) {
		this.blockStarts = blockStarts;
	}

	public String getBlockStarts() {
		return blockStarts;
	}

	
	
	
}
