package ch.epfl.bbcf.parser;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ch.epfl.bbcf.exception.ParsingException;
import ch.epfl.bbcf.feature.Feature;
import ch.epfl.bbcf.feature.Track;


public abstract class Parser {

	public enum Processing {SEQUENCIAL,TOTAL};

	protected Map<Track,List<Feature>> features;
	protected Processing processingType;

	/**
	 * initialize a new parser with the type of processing
	 * and it's handler
	 * @param type
	 * @param handler
	 */
	protected Parser(Processing type){
		switch(type){
		case TOTAL :
			features = new HashMap<Track,List<Feature>>();
			break;
		}
		processingType = type;
		lineNb=0;
	}

	/**
	 * get the features of the file :
	 * key = the track; value = list of features
	 * this will be null if you have choosen the SEQUENCIAL processing type
	 * @return a Map
	 */
	public Map<Track, List<Feature>> getFeatures() {
		return features;
	}
	
	/**
	 * the current line being parsed
	 */
	protected int lineNb;

	/**
	 * method called when parsing the file
	 * @param line - the current line
	 */
	protected abstract void processLine(String line,Handler handler) throws ParsingException;

	/**
	 * method to call to parse the file
	 * @param file - the file
	 * @throws IOException
	 * @throws ParsingException 
	 */
	public void parse(File file,Handler handler) throws IOException, ParsingException{
		FileReader fr = null;
		BufferedReader br = null;
		fr = new FileReader(file);
		br = new BufferedReader(fr);
		String line;
		handler.start();
		while ((line = br.readLine()) != null){
			lineNb++;
			processLine(line,handler);
		}
		handler.end();
	}
	/**
	 * method called when a new feature is parsed
	 * @param feature
	 */
	protected void newFeature(Handler handler,Track track,Feature feature){
		switch(processingType){
		case TOTAL :
			List<Feature> feats = features.get(track);
			feats.add(feature);
			features.put(track, feats);
			break;
		case SEQUENCIAL :
			handler.newFeature(feature);
			break;
		}
	}
	/**
	 * method called when a newTrack is parsed
	 * @param track
	 */
	protected void newTrack(Handler handler,Track track){
		switch(processingType){
		case TOTAL:
			features.put(track, new ArrayList<Feature>());
			break;
		case SEQUENCIAL:
			handler.newTrack(track);
			break;
		}
	}
	protected int getInt(String str) throws ParsingException {
		try {
			return Integer.parseInt(str);
		} catch(NumberFormatException nfe){
			throw new ParsingException(nfe,"NumberFormatException",lineNb);
		}
	}
	protected float getScore(String s) throws ParsingException {
		try{
			return Float.parseFloat(s);
		} catch (NumberFormatException e){
			throw new ParsingException(e,"NumberFormatException",lineNb);
		}
	}
	protected int getStrand(String strand) throws ParsingException {
		if(strand.equalsIgnoreCase("+") || strand.equalsIgnoreCase("1")){
			return 1;
		} else if(strand.equalsIgnoreCase("-") || strand.equalsIgnoreCase("-1")){
			return -1;
		} else {
			throw new ParsingException("Strand "+strand+" not recognized",lineNb);
		}
	}
	

}
