package ch.epfl.bbcf.parser;


import ch.epfl.bbcf.feature.Feature;
import ch.epfl.bbcf.feature.Track;



public interface Handler {

	/**
	 * called when a new feature is encountered in the document
	 * @param feature
	 * @throws JSONException 
	 */
	public abstract void newFeature(Feature feature);

	/**
	 * called when a new track is encountered in the document
	 * @param track
	 */
	public abstract void newTrack(Track track);
	
	/**
	 * method called at the start of the file
	 */
	public abstract void start();
	/**
	 * method called at EOF
	 */
	public abstract void end();

	
}
