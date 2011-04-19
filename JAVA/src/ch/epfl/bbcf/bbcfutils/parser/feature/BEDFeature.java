package ch.epfl.bbcf.bbcfutils.parser.feature;

public class BEDFeature extends QualitativeFeature{

	private Integer thickStart,thickEnd;
	private String itemRgb,blockCount,blockSize,blockStarts;
	
	public BEDFeature(String chromosome2, Integer start2, Integer end2, String name2,
			Integer strand2, Float score2, Integer thickStart, Integer thickEnd, 
			String itemRgb, String blockCount, String blockSizes, String blockStarts) {
		super(chromosome2, start2, end2, score2, strand2, name2,"");
		this.attributes = buildAttributesfromBEDFeature(
				thickEnd, thickEnd, itemRgb, 
				blockCount, blockSizes, blockStarts);
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

	public void setThickStart(Integer thickStart) {
		this.thickStart = thickStart;
	}

	public Integer getThickStart() {
		return thickStart;
	}

	public void setThickEnd(Integer thickEnd) {
		this.thickEnd = thickEnd;
	}

	public Integer getThickEnd() {
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
