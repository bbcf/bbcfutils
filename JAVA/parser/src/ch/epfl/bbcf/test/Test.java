package ch.epfl.bbcf.test;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import ch.epfl.bbcf.exception.ParsingException;
import ch.epfl.bbcf.feature.Feature;
import ch.epfl.bbcf.feature.Track;
import ch.epfl.bbcf.feature.WIGFeature;
import ch.epfl.bbcf.parser.Parser;
import ch.epfl.bbcf.parser.WIGParser;
import ch.epfl.bbcf.parser.Parser.Processing;

public class Test{


	/**
	 * this handler will give you access to the file
	 * @author Yohan Jarosz
	 *
	 */
	public static class MyParsingHandler implements ch.epfl.bbcf.parser.Handler{
		@Override
		public void start() {
			System.out.println("start of parsing");
		}
		@Override
		public void end() {
			System.out.println("end of parsing");
		}
		@Override
		public void newFeature(Feature feature) {
			System.out.println("new feature");
			WIGFeature wigFeat = (WIGFeature) feature;
			System.out.println(wigFeat.detail());
		}
		@Override
		public void newTrack(Track track) {
			System.out.println("new track");
			System.out.println(track.detail());
		}


	}
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		System.out.println("STARTING OF SEQUENCIAL PARSING");
		/**
		 * in order to parse hudge files,
		 * nothing is kept in memory,
		 * so each time a new feature is encountered, 
		 * newFeature(Feature feature) is called
		 * (newTrack(Track track) when a new track ... )
		 */
		long start = System.currentTimeMillis();
		Parser parser = new WIGParser(Processing.SEQUENCIAL);
		MyParsingHandler handler = new Test.MyParsingHandler();
		try {
			parser.parse(new File(
			"/path/to/wiggle/file") ,handler);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (ParsingException e) {
			e.printStackTrace();
		}
		long end = System.currentTimeMillis();
		long time = (end-start);
		System.out.println(" - "+time);		

		
		
		System.out.println("STARTING OF TOTAL PARSING");
		/**
		 * everything is parsed an then you have access to the features
		 * via  Map<Track, List<Feature>> features = parser.getFeatures();
		 */
		start = System.currentTimeMillis();
		parser = new WIGParser(Processing.TOTAL);
		try {
			parser.parse(new File(
			"/path/to/wiggle/file") ,handler);
		} catch (IOException e) {
			e.printStackTrace();
		} catch (ParsingException e) {
			e.printStackTrace();
		}

		end = System.currentTimeMillis();
		time = (end-start);
		System.out.println(" - "+time);	
		
		Map<Track, List<Feature>> feats = parser.getFeatures();
		for(Map.Entry<Track, List<Feature>> entry : feats.entrySet()){
			System.out.println(entry.getKey().detail());
			for(Feature feat : entry.getValue()){
				WIGFeature wigFeat = (WIGFeature) feat;
				System.out.println(wigFeat.detail());
			}
			System.out.println();
		}

	}
}
