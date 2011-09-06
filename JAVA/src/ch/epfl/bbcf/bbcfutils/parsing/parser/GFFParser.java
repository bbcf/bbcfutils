package ch.epfl.bbcf.bbcfutils.parsing.parser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.biojava3.genome.parsers.gff.Location;


import ch.epfl.bbcf.bbcfutils.exception.ParsingException;
import ch.epfl.bbcf.bbcfutils.parsing.feature.GFFFeature;
import ch.epfl.bbcf.bbcfutils.parsing.parser.Parser.Processing;


/**
 * this class parse GFF files - 
 * it use the GFFParser from biojava -
 * if you want more control over the file, you should use biojava parser instead
 * @author Yohan Jarosz
 *
 */
public class GFFParser extends Parser{

	private static String[] hierarchy = {"id","gene_id"};
	private static String[] display = {"name","gene_name"};
	private int lineNb;
	private Handler handler;
	private List<String> types;

	public GFFParser(Processing type) {
		super(type);
		setTypes(new ArrayList<String>());
		lineNb=0;
	}

	public void parse(File gffFile,Handler handler) throws IOException, ParsingException{
		this.setHandler(handler);
		BufferedReader br = new BufferedReader(new FileReader(gffFile));
		String s;
		handler.start();
		for (s = br.readLine(); null != s; s = br.readLine()) {
			s = s.trim();
			if (s.length() > 0) {
				if (s.charAt(0) == '#') {
					//ignore comment lines
					if(s.startsWith("##fasta"))
						break;
				} else {
					GFFFeature f = parseLine(s);
					if (f != null) {
						handler.newFeature(f);
					}
				}
			}

		}
		br.close();
		handler.end();
	}

	/**
	 * create Feature from line of GFF file
	 */
	private GFFFeature parseLine(String s) {
		lineNb++;
		//FIXME update to use regex split on tabs
		//FIXME better errors on parse failures
		int start = 0;
		int end = 0;

		start = end;
		end = s.indexOf('\t', start);
		String seqname = s.substring(start, end).trim();

		start = end + 1;
		end = s.indexOf('\t', start);
		//TODO use source to put right start & end
		String source = s.substring(start, end).trim();

		start = end + 1;
		end = s.indexOf('\t', start);
		String type = s.substring(start, end);

		start = end + 1;
		end = s.indexOf('\t', start);
		String locStart = s.substring(start, end);

		start = end + 1;
		end = s.indexOf('\t', start);
		String locEnd = s.substring(start, end);

		float score;
		start = end + 1;
		end = s.indexOf('\t', start);
		try {
			score = Float.parseFloat(s.substring(start, end));
		} catch (Exception e) {
			score = 0f;
		}

		start = end + 1;
		end = s.indexOf('\t', start);
		char strand = s.charAt(end - 1);
		String strand_str =  s.substring(start, end);
		//added by scooter willis to deal with glimmer predictions that
		//have the start after the end but is a negative strand
		int locationStart = Integer.parseInt(locStart);
		int locationEnd = Integer.parseInt(locEnd);
		if(locationStart > locationEnd){
			int temp = locationStart;
			locationStart = locationEnd;
			locationEnd = temp;

		}
		int strand_int =0;
		if(strand_str.equalsIgnoreCase("-")){
			strand_int=-1;
		} else if(strand_str.equalsIgnoreCase("+")){
			strand_int=+1;
		} else {
			try {
				strand_int = Integer.parseInt(strand_str);
			} catch (Exception e) {
				strand_int = 0;
			}
		}
		
		
		
		
		
		Location location = Location.fromBio(locationStart, locationEnd, strand);
		assert (strand == '-') == location.isNegative();
		int frame;
		//FIXME 
		start = end + 1;
		end = s.indexOf('\t', start);
		String f = s.substring(start, end);
		

		//grab everything until end of line (or # comment)
		start = end + 1;
		end = s.indexOf('#', start);
		String attributes = null;
		if (end < 0) {
			attributes = new String(s.substring(start));
		} else {
			attributes = new String(s.substring(start, end));
		}
		if(!types.contains(type)){
			types.add(type);
		}
		String displayName="";
		String hierarchyId="";
		Map<String,String> attMap = takeAttributesMap(attributes);
		for(String hier:hierarchy){
			if(attMap.containsKey(hier)){
				hierarchyId=attMap.get(hier);
				break;
			}
		}
		for(String dis:display){
			if(attMap.containsKey(dis)){
				displayName=attMap.get(dis);
				break;
			}
		}
		return new GFFFeature(seqname,source,locationStart,locationEnd,
				displayName,hierarchyId,score, strand_int, type, attributes);

	}

	private static Map<String, String> takeAttributesMap(String attributes) {
		Map<String, String> attMap = new HashMap<String, String>();
		String[]atts = attributes.split(";");
		for(String att:atts){
			att = att.trim();
			String[] keyvalue = att.split("\\s");
			if(keyvalue.length>1){
				attMap.put(keyvalue[0],keyvalue[1]);
			}
		}
		return attMap;
	}

	@Override
	protected void processLine(String line, Handler handler)
	throws ParsingException {
		//not used
	}

	public void setHandler(Handler handler) {
		this.handler = handler;
	}

	public Handler getHandler() {
		return handler;
	}

	public void setTypes(List<String> types) {
		this.types = types;
	}

	public List<String> getTypes() {
		return types;
	}

	
	

}
