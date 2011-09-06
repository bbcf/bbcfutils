package ch.epfl.bbcf.bbcfutils.parsing.feature;



public class WIGTrack extends Track{


	public static final int WIG_TRACK_TYPE = 1;
	public static final int BEDGRAPH_TRACK_TYPE = 2;
	private int trackType;
	
	
	
	
	
	
	
	public WIGTrack(){
		super();
	}
	
	public void setTrackType(int trackType) {
		this.trackType = trackType;
	}
	/**
	 * get the track type : (WIG or Bedgraph)
	 * @return 1 for wig, 2 for bedgraph
	 */
	public int getTrackType() {
		return trackType;
	}
	

	

	
	
}
