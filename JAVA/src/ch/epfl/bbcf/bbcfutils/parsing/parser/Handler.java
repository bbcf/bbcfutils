package ch.epfl.bbcf.bbcfutils.parsing.parser;

import java.sql.SQLException;

import ch.epfl.bbcf.bbcfutils.exception.ParsingException;
import ch.epfl.bbcf.bbcfutils.parsing.feature.Feature;
import ch.epfl.bbcf.bbcfutils.parsing.feature.Track;





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
	 * @throws ParsingException 
	 */
	public abstract void newTrack(Track track) throws ParsingException;
	
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
