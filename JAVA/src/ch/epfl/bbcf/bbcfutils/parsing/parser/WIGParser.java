package ch.epfl.bbcf.bbcfutils.parsing.parser;


import java.util.regex.Matcher;
import java.util.regex.Pattern;

import ch.epfl.bbcf.bbcfutils.exception.ParsingException;
import ch.epfl.bbcf.bbcfutils.parsing.feature.WIGFeature;
import ch.epfl.bbcf.bbcfutils.parsing.feature.WIGTrack;







public class WIGParser extends Parser{

	/**
	 * the current track
	 */
	private WIGTrack cur_track;
	/**
	 * the current feature
	 */
	private WIGFeature cur_feature;
	/**
	 * the previous feature
	 */
	private WIGFeature prev_feature;

	/**
	 * if the track parameters are finished to read
	 */
	private boolean trackParametersRead;

	/**
	 * the lines beginning by 'browser'
	 */
	private String browserLine;

	/**
	 * a boolean to know if there is a new track
	 */
	private boolean newTrack;

	/**
	 * pattern matching the track's attributes
	 */
	private static final Pattern trackAttributesPattern = 
		Pattern.compile("(\\w+=\"\\w+ \\w+\")|(\\w+=\"\\w+\")|(\\w+=\\w+)|(\\w+=.+)");

	public WIGParser(Processing type) {
		super(type);
		trackParametersRead = false;
		this.browserLine="";
	}


	@Override
	protected void processLine(String line, Handler handler)
	throws ParsingException {
		if(line.startsWith("#")){
			//skip comments
		}else if(line.startsWith("track type=")){
			newTrack=true;
			//new track
			trackParametersRead = false;
			cur_track = new WIGTrack();
			cur_feature=null;
			prev_feature=null;
			handleTrackType(line);
		} else if(line.startsWith("browser")){
			addToBrowserLine(line);
		} else if(line.startsWith("fixedStep")||line.startsWith("variableStep")){
			if(newTrack){
				newTrack=false;
				newTrack(handler, cur_track);
			}
			//new feature
			prev_feature = null;
			cur_feature = new WIGFeature();
			trackParametersRead = true;
			handleFormatLine(line);
		} else {
			if(trackParametersRead){
				switch(cur_track.getTrackType()){
				case WIGTrack.WIG_TRACK_TYPE:handleWIGLine(line,handler);break;
				case WIGTrack.BEDGRAPH_TRACK_TYPE:handlerBEDLine(line,handler);break;
				}
			} else {
				String[] tab = line.split("\\s");
				if(tab.length>1){
					if(null!=Integer.valueOf(tab[2]));
					trackParametersRead = true;
					processLine(line, handler);
				} else {
					handleTrackType(line);
				}

			}

		}

	}

	/**
	 * add the parameter to the browser line
	 * @param line
	 */
	private void addToBrowserLine(String line) {
		this.browserLine+=line.substring(8, line.length())+" ";
	}
	/**
	 * return the 'browser' line of a wig file
	 * @return the first line of the wig file
	 */
	public String getBrowserLine() {
		return browserLine;
	}
	/**
	 * handle the line in BED format
	 * @param line
	 * @param handler 
	 * @throws ParsingException 
	 */
	private void handlerBEDLine(String line, Handler handler) throws ParsingException {
		String[] chr_start_end_score = line.split("\\s");
		Integer start,end;
		Float score;
		String chr;
		switch(chr_start_end_score.length){
		case 4: 
			chr = chr_start_end_score[0];
			start = getInt(chr_start_end_score[1]);
			end = getInt(chr_start_end_score[2]);
			score = getScore(chr_start_end_score[3]);
			cur_feature = new WIGFeature(chr,start,end,score);
			break;
		default: throw new ParsingException("the entry don't have required number of fields " +
				"(4 required : chromosome,start,end,score separated by spaces or tabs): ", lineNb);
		}
		newFeature(handler, cur_track, cur_feature);
	}
	/**
	 * handle the line in wiggle format
	 * @param line
	 * @param handler 
	 * @throws ParsingException 
	 */
	private void handleWIGLine(String line, Handler handler) throws ParsingException {
		switch(cur_feature.getFormat()){
		case WIGFeature.FIXED_STEP:
			cur_feature.setScore(getScore(line));
			if(null!=prev_feature){
				cur_feature.setStart(prev_feature.getStart()+prev_feature.getStep());
				cur_feature.setEnd(prev_feature.getEnd()+prev_feature.getStep());
				cur_feature.setChromosome(prev_feature.getChromosome());
			} else {
				cur_feature.setEnd(cur_feature.getStart()+cur_feature.getSpan());
			}
			break;
		case WIGFeature.VARIABLE_STEP:
			String[] pos_value = line.split("\\s");
			cur_feature.setStart(getInt(pos_value[0]));
			cur_feature.setScore(getScore(pos_value[1]));
			cur_feature.setEnd(cur_feature.getStart()+cur_feature.getSpan());
			if(null!=prev_feature){
				cur_feature.setChromosome(prev_feature.getChromosome());
			}
			break;
		}
		newFeature(handler, cur_track,new WIGFeature(cur_feature));
		prev_feature = cur_feature;
	}
	/**
	 * handle the line beginning with fixedstep or 
	 * variablestep
	 * @param line
	 * @throws ParsingException 
	 */
	private void handleFormatLine(String line) throws ParsingException {
		String[] params = line.split("\\s");
		for (int i = 0; i < params.length; i++) {
			String[] str = params[i].split("=");
			if(str[0].startsWith("fixedStep")){
				cur_feature.setFormat(WIGFeature.FIXED_STEP);
			} else if(str[0].startsWith("variableStep")){
				cur_feature.setFormat(WIGFeature.VARIABLE_STEP);
			} else if (str[0].equalsIgnoreCase("chrom")) {
				cur_feature.setChromosome(str[1]);
			} else if (str[0].equalsIgnoreCase("span")) {
				try{
					cur_feature.setSpan(Integer.valueOf(str[1]));
				}catch (NumberFormatException e){
					throw new ParsingException("NumberFormatException - the span must be an integer - ",lineNb);
				}
			} else if (str[0].equalsIgnoreCase("step")) {
				try {
					cur_feature.setStep(Integer.valueOf(str[1]));
				}catch (NumberFormatException e){
					throw new ParsingException("NumberFormatException - the step must be an integer - ",lineNb);
				}
			} else if (str[0].equalsIgnoreCase("start")) {
				try {
					cur_feature.setStart(Integer.valueOf(str[1]));
				}catch (NumberFormatException e){
					throw new ParsingException("NumberFormatException - the start must be an integer - ",lineNb);
				}
			}
		}
		if(null==cur_feature.getSpan()){
			cur_feature.setSpan(1);
		}
	}

	/**
	 * handle the line REQUIRED in a wiggle file
	 * to determine the track type
	 * @param line
	 * @throws ParsingException 
	 */
	private void handleTrackType(String line) throws ParsingException {
		String[] tmp = line.split("=");
		if (tmp[1].startsWith("bedGraph")) {
			cur_track.setTrackType(WIGTrack.BEDGRAPH_TRACK_TYPE);
		} else if (tmp[1].startsWith("wiggle_0")) {
			cur_track.setTrackType(WIGTrack.WIG_TRACK_TYPE);
		}
		Matcher m = trackAttributesPattern.matcher(line);
		while(m.find()){
			String[]tab = m.group().split("=");
			cur_track.addAttribute(tab[0],tab[1]);
		}
	}

}

