package ch.epfl.bbcf.bbcfutils.conversion.sqlite;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import ch.epfl.bbcf.bbcfutils.access.genrep.GenrepWrapper;
import ch.epfl.bbcf.bbcfutils.access.genrep.MethodNotFoundException;
import ch.epfl.bbcf.bbcfutils.access.genrep.json_pojo.Assembly;
import ch.epfl.bbcf.bbcfutils.access.genrep.json_pojo.Chromosome;
import ch.epfl.bbcf.bbcfutils.exception.ExtensionNotRecognisedException;
import ch.epfl.bbcf.bbcfutils.parser.BEDParser;
import ch.epfl.bbcf.bbcfutils.parser.GFFParser;
import ch.epfl.bbcf.bbcfutils.parser.Handler;
import ch.epfl.bbcf.bbcfutils.parser.Parser;
import ch.epfl.bbcf.bbcfutils.parser.WIGParser;
import ch.epfl.bbcf.bbcfutils.parser.Parser.Processing;
import ch.epfl.bbcf.bbcfutils.parser.exception.ParsingException;
import ch.epfl.bbcf.bbcfutils.parser.feature.ExtendedQualitativeFeature;
import ch.epfl.bbcf.bbcfutils.parser.feature.Feature;
import ch.epfl.bbcf.bbcfutils.parser.feature.QualitativeFeature;
import ch.epfl.bbcf.bbcfutils.parser.feature.QuantitativeFeature;
import ch.epfl.bbcf.bbcfutils.parser.feature.Track;
import ch.epfl.bbcf.bbcfutils.sqlite.SQLiteAccess;
import ch.epfl.bbcf.bbcfutils.sqlite.SQLiteConstruct;

public class ConvertToSQLite {

	public enum Extension {WIG,BEDGRAPH,GFF,BED,BAM}
	private String inputPath;

	private Parser parser;
	private Handler handler;

	private Extension extension;
	private SQLiteConstruct construct;
	//if an nr assembly is provided
	//we have to check the names
	private int nrAssemblyId;
	private List<String> chromosomes;
	private Map<String,String> altsNames;
	private String previousUnmapped;

	private Assembly assembly;

	private String outputPath;

	//hash map to fill types attributes for an extended track
	private List<String> types;

	public ConvertToSQLite(String inputPath,Extension extension) throws ExtensionNotRecognisedException{
		this.inputPath = inputPath;
		this.parser = takeParser(extension);
		this.handler = takeHandler();
		this.extension=extension;
		this.nrAssemblyId=-1;
	}
	public ConvertToSQLite(String inputPath,Extension extension,int nrAssemblyId) throws ExtensionNotRecognisedException, ParsingException{
		System.out.println("cv to SQlite");
		this.inputPath = inputPath;
		this.parser = takeParser(extension);
		this.handler = takeHandler();
		this.extension=extension;
		this.nrAssemblyId=nrAssemblyId;
		try {
			this.chromosomes = takeChromosomes(nrAssemblyId);
		} catch (MethodNotFoundException e) {
			throw new ParsingException(e);
		} catch (IOException e) {
			throw new ParsingException(e);
		}
		this.altsNames=new HashMap<String, String>();
		this.previousUnmapped="";
		this.assembly = takeAssembly(nrAssemblyId);
	}

	private Assembly takeAssembly(int nrAssemblyId2) {
		try {
			return GenrepWrapper.getAssemblyFromNrAssemblyId(nrAssemblyId);
		} catch (MethodNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}
	/**
	 * fetch the list of chromosomes from Genrep
	 * @param nrAssemblyId2 the nr assembly identifier
	 * @return a list of chromosomes names
	 * @throws IOException 
	 * @throws MethodNotFoundException 
	 */
	private List<String> takeChromosomes(int nrAssemblyId2) throws MethodNotFoundException, IOException {
		Assembly assembly;
		System.out.println("take chromosomes");
		assembly = GenrepWrapper.getAssemblyFromNrAssemblyId(nrAssemblyId2);
		List<Chromosome> chromosomes = assembly.getChromosomes();
		List<String> chrNames = new ArrayList<String>();
		for(Chromosome chromosome : chromosomes){
			System.out.println("chromosome : "+chromosome.getChr_name());
			chrNames.add(chromosome.getChr_name());
		}
		return chrNames;
	}

	/**
	 * initialize the parsing handler
	 * @return a ParsingHandler
	 */
	private Handler takeHandler() {
		return new ParsingHandler();
	}


	/**
	 * get the right parser for the 
	 * right extension
	 * @param extension - the extension
	 * @return the parser
	 * @throws ExtensionNotRecognisedException 
	 */
	private static Parser takeParser(Extension extension) throws ExtensionNotRecognisedException {
		Parser p = null;
		switch(extension){
		case WIG:case BEDGRAPH:
			p = new WIGParser(Processing.SEQUENCIAL);
			break;
		case BED:
			p = new BEDParser(Processing.SEQUENCIAL);
			break;
		case GFF:
			p = new GFFParser(Processing.SEQUENCIAL);
			break;
			//		case BAM:
			//			p = new BAMParser(Processing.SEQUENCIAL);
			//			break;
		default :
			throw new ExtensionNotRecognisedException(extension);
		}
		return p;
	}


	public void setInputPath(String inputPath) {
		this.inputPath = inputPath;
	}
	public String getInputPath() {
		return inputPath;
	}


	/**
	 * launch the conversion
	 * @param outputPath - where the output should go
	 * @return true if successful
	 * @throws IOException
	 * @throws ParsingException
	 * @throws InstantiationException
	 * @throws IllegalAccessException
	 * @throws ClassNotFoundException
	 * @throws SQLException
	 */
	public boolean convert(String outputPath,String type) throws IOException, ParsingException{
		this.outputPath = outputPath;
		File input = new File(inputPath);
		try {
			this.construct = SQLiteConstruct.getConnectionWithDatabase(outputPath);
			this.construct.createNewDatabase(type);
			this.construct.commit();
		} catch (InstantiationException e) {
			throw new ParsingException(e);
		} catch (IllegalAccessException e) {
			throw new ParsingException(e);
		} catch (ClassNotFoundException e) {
			throw new ParsingException(e);
		} catch (SQLException e) {
			throw new ParsingException(e);
		}
		if(!input.exists()){
			throw new FileNotFoundException(inputPath);
		}
		//initialize the map for extended type
		switch(extension){
		case GFF :
			this.types=new ArrayList<String>();
		}
		this.parser.parse(input, handler);
		return true;
	}











	/**
	 * Class which handle the parsing of the file
	 * @author Yohan Jarosz
	 *
	 */
	private class ParsingHandler implements Handler{
		@Override
		public void newFeature(Feature feature) throws ParsingException {
			String chromosome = feature.getChromosome();
			//change the chromosome name if a 
			//nr assembly id is provided
			try {
				chromosome=guessChromosome(chromosome,nrAssemblyId);
			} catch (MethodNotFoundException e1) {
				throw new ParsingException(e1);
			} catch (IOException e1) {
				throw new ParsingException(e1);
			}
			if(null==chromosome){
				return;
			}
			try {
				switch(extension){
				case GFF : 
					ExtendedQualitativeFeature feat = (ExtendedQualitativeFeature)feature;
					if(!construct.isCromosomeCreated(chromosome)){
						construct.newChromosome_qual_extended(chromosome);
					}
					int index = types.indexOf(feat.getType());
					if(-1==index){
						types.add(feat.getType());
						index=types.size();
					}
					construct.writeValues_qual_extended(
							chromosome, feat.getStart(), 
							feat.getEnd(),feat.getScore(), 
							feat.getName(),feat.getStrand(),
							feat.getAttributes(),
							index,
							feat.getIdentifier());
					break;
				case WIG:
					QuantitativeFeature feat1 = (QuantitativeFeature)feature;
					if(!construct.isCromosomeCreated(chromosome)){
						construct.newChromosome_quant(chromosome);
					}

					construct.writeValues_quant(chromosome, feat1.getStart(), 
							feat1.getEnd(), feat1.getScore());
					break;
				case BED: case BAM:
					QualitativeFeature feat2 = (QualitativeFeature)feature;
					if(!construct.isCromosomeCreated(chromosome)){
						construct.newChromosome_qual(chromosome);
					}
					construct.writeValues_qual(
							chromosome, feat2.getStart(), feat2.getEnd(), 
							feat2.getScore(),
							feat2.getName(),
							feat2.getStrand(),
							feat2.getAttributes());
					break;
				}
			} catch (SQLException e) {
				throw new ParsingException(e);
			}
		}

		private String guessChromosome(String current, Integer nrAssemblyId) throws MethodNotFoundException, IOException {
			if(nrAssemblyId!=null){
				//if current == previous unfinded, return null
				if(current.equalsIgnoreCase(previousUnmapped)){
					return null;
				}
				//if current is already a good one, return it
				if(chromosomes.contains(current)){
					return current;
				}
				//if mapped alredy found, return the mapped
				if(altsNames.containsKey(current)){
					return altsNames.get(current);
				}
				//try to find a mapping
				Chromosome newChr = GenrepWrapper.guessChromosome(current, assembly.getId());

				if(null==newChr){//no mapping found
					previousUnmapped=current;
					return null;
				} else {//mapping founded
					String mapping =newChr.getChr_name();
					altsNames.put(current, mapping);
					return mapping;
				}
			}
			return current;
		}
		@Override
		public void newTrack(Track track) {
			System.err.println("Operation not supported");
		}

		@Override
		public void start() throws ParsingException {
			System.out.println("start");
			try {
				switch(extension){
				case WIG:
					construct.createNewDatabase("quantitative");
					break;
				case BED:case BAM:case GFF:
					construct.createNewDatabase("qualitative");
					break;
				}
			} catch (SQLException e) {
				throw new ParsingException(e);
			}
		}

		@Override
		public void end() throws ParsingException {
			System.out.println("end");
			try {
				construct.commit();
				List<String> chrNames = construct.getChromosomesNames();
				if(chrNames.isEmpty()){
					throw new ParsingException("no chromosomes found in building the database");
				}
				Map<String,Integer> map = new HashMap<String, Integer>();
				for(String chr : chrNames){
					SQLiteAccess access = SQLiteAccess.getConnectionWithDatabase(outputPath);
					int length = access.getMaxEndForChromosome(chr);
					access.close();
					if(length!=0){
						map.put(chr,length);
					}
				}
				switch(extension){
				case WIG:
					construct.finalizeDatabase(map, false, false, true);
					break;
				case BED:case BAM:
					construct.finalizeDatabase(map, true, false, true);
					break;
				case GFF:
					construct.finalizeExtended_qual(types);
					construct.finalizeDatabase(map, true, false, true);
					break;
				}
				construct.commit();
				construct.close();
			} catch (SQLException e) {
				throw new ParsingException(e);
			} catch (InstantiationException e) {
				throw new ParsingException(e);
			} catch (IllegalAccessException e) {
				throw new ParsingException(e);
			} catch (ClassNotFoundException e) {
				throw new ParsingException(e);
			}
			System.out.println("conversion to SQLite done ");
		}

	}



	public void setNrAssemblyId(int nrAssemblyId) {
		this.nrAssemblyId = nrAssemblyId;
	}


	public int getNrAssemblyId() {
		return nrAssemblyId;
	}
	public void setChromosomes(List<String> chromosomes) {
		this.chromosomes = chromosomes;
	}
	public List<String> getChromosomes() {
		return chromosomes;
	}

	public static void main(String[] args) throws MethodNotFoundException, IOException{
		//		try {
		//			ConvertToSQLite c = new ConvertToSQLite("/Users/jarosz/Desktop/toto.gtf",Extension.GFF,70);
		//			c.convert("/Users/jarosz/Desktop/A4.db","qualitative");
		//		} catch (IOException e) {
		//			// TODO Auto-generated catch block
		//			e.printStackTrace();
		//		} catch (ParsingException e) {
		//			// TODO Auto-generated catch block
		//			e.printStackTrace();
		//		} catch (ExtensionNotRecognisedException e) {
		//			// TODO Auto-generated catch block
		//			e.printStackTrace();
		//		}



		String previousUnmapped ="";
		String chromosome = "MT";
		Integer nrAssemblyId = 70;
		Map<String, String> altsNames=new HashMap<String, String>();
		List<String> chromosomes = new ArrayList<String>();
		Assembly assembly = GenrepWrapper.getAssemblyFromNrAssemblyId(nrAssemblyId);
		List<Chromosome> chromosomes_ = assembly.getChromosomes();
		for(Chromosome chr : chromosomes_){
			chromosomes.add(chr.getChr_name());
		}
		System.out.println(chromosomes);

		if(nrAssemblyId!=null){
			if(!chromosome.equalsIgnoreCase(previousUnmapped)){
				if(!chromosomes.contains(chromosome)){
					if(altsNames.containsKey(chromosome)){
						chromosome=altsNames.get(chromosome);
					} else {
						System.out.println("gzess chr");
						Chromosome newChr = GenrepWrapper.guessChromosome(chromosome, assembly.getId());
						System.out.println(" "+newChr);
						if(null==newChr){
							previousUnmapped=chromosome;
							System.out.println("NULL1");
						} else {
							String tmp =newChr.getChr_name();
							System.out.println("J�IJ�OIJ "+tmp);
							if(!chromosomes.contains(tmp)){
								altsNames.put(chromosome, tmp);
								chromosome=tmp;
							} else {

								System.out.println("NULL2");
							}
						}
					}
				}
			} else {
				System.out.println("NULL3");
			}
		}
		System.out.println(chromosome);












	}
}
