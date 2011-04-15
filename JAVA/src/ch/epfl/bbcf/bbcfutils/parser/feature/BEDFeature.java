package ch.epfl.bbcf.bbcfutils.parser.feature;

public class BEDFeature extends QualitativeFeature{

	private int thickStart,thickEnd;
	private String itemRgb,blockCount,blockSize,blockStarts;
	
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
	
	

	@Override
	public String detail(){
		return super.detail()+"\n"+
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
