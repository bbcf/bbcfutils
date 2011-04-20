package ch.epfl.bbcf.bbcfutils.parser;

import ch.epfl.bbcf.bbcfutils.parser.exception.ParsingException;
import ch.epfl.bbcf.bbcfutils.parser.feature.Track;
import ch.epfl.bbcf.bbcfutils.parser.feature.Feature;





public interface Handler {

	/**
	 * called when a new feature is encountered in the document
	 * @param feature
	 * @throws ParsingException 
	 * @throws JSONException 
	 */
	public abstract void newFeature(Feature feature) throws ParsingException;

	/**
	 * called when a new track is encountered in the document
	 * @param track
	 */
	public abstract void newTrack(Track track);
	
	/**
	 * method called at the start of the file
	 * @throws ParsingException 
	 */
	public abstract void start() throws ParsingException;
	/**
	 * method called at EOF
	 * @throws ParsingException 
	 */
	public abstract void end() throws ParsingException;

	
}
