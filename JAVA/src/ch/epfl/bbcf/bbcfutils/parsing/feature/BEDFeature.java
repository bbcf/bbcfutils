package ch.epfl.bbcf.bbcfutils.parsing.feature;

public class BEDFeature extends Feature{

	private String chromosome,name;
	

	private Integer start,end,strand,thickStart,thickEnd;
	private String itemRgb,blockCount,blockSize,blockStarts;
	private Float score;
	
	public BEDFeature(String chromosome2, Integer start2, Integer end2, String name2,
			Integer strand2, Float score2, Integer thickStart, Integer thickEnd, 
			String itemRgb, String blockCount, String blockSizes, String blockStarts) {
		this.chromosome=chromosome2;
		this.start=start2;
		this.end = end2;
		this.name=name2;
		this.strand=strand2;
		this.score=score2;
		this.thickStart=thickStart;
		this.thickEnd=thickEnd;
		this.itemRgb=itemRgb;
		this.blockCount=blockCount;
		this.blockSize=blockSizes;
		this.blockStarts=blockStarts;
	}
	
	@Override
	public String display(){
		return ds(this.chromosome)+SEP+ds(this.start)
		+SEP+ds(this.end)+SEP+ds(this.name)+SEP+ds(this.score)
		+SEP+dss(this.strand)+SEP+ds(this.thickStart)+
		SEP+ds(this.thickEnd)+SEP+ds(this.itemRgb)+
		SEP+ds(this.blockCount)+SEP+ds(this.blockSize)+
		SEP+ds(this.blockStarts);
		
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

	
	public String getChromosome() {
		return chromosome;
	}

	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}

	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}

	public Integer getStart() {
		return start;
	}

	public void setStart(Integer start) {
		this.start = start;
	}

	public Integer getEnd() {
		return end;
	}

	public void setEnd(Integer end) {
		this.end = end;
	}

	public Integer getStrand() {
		return strand;
	}

	public void setStrand(Integer strand) {
		this.strand = strand;
	}

	public Float getScore() {
		return score;
	}

	public void setScore(Float score) {
		this.score = score;
	}
	
}
