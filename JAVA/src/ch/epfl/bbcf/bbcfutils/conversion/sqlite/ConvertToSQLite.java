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

	public ConvertToSQLite(String inputPath,Extension extension) throws ExtensionNotRecognisedException{
		this.inputPath = inputPath;
		this.parser = takeParser(extension);
		this.handler = takeHandler();
		this.extension=extension;
		this.nrAssemblyId=-1;
	}
	public ConvertToSQLite(String inputPath,Extension extension,int nrAssemblyId) throws ExtensionNotRecognisedException, ParsingException{
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
		assembly = GenrepWrapper.getAssemblyFromNrAssemblyId(nrAssemblyId2);
		List<Chromosome> chromosomes = assembly.getChromosomes();
		List<String> chrNames = new ArrayList<String>();
		for(Chromosome chromosome : chromosomes){
			chrNames.add(chromosome.getName());
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
					construct.writeValues_qual_extended(
							chromosome, feat.getStart(), 
							feat.getEnd(),feat.getScore(), 
							feat.getName(),feat.getStrand(),
							feat.getAttributes(),
							feat.getType(),
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

		private String guessChromosome(String chromosome, Integer nrAssemblyId) throws MethodNotFoundException, IOException {
			if(nrAssemblyId!=null){
				if(!chromosome.equalsIgnoreCase(previousUnmapped)){
					if(!chromosomes.contains(chromosome)){
						if(altsNames.containsKey(chromosome)){
							chromosome=altsNames.get(chromosome);
						} else {

							Chromosome newChr = GenrepWrapper.guessChromosome(chromosome, assembly.getId());
							if(null==newChr){
								previousUnmapped=chromosome;
								return null;
							} else {
								String tmp =newChr.getName();
								altsNames.put(chromosome, tmp);
								chromosome=tmp;
							}
						}
					}
				} else {
					return null;
				}
			}
			return chromosome;
		}

		@Override
		public void newTrack(Track track) {
			System.err.println("Operation not supported");
		}

		@Override
		public void start() throws ParsingException {
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
			try {
				construct.commit();
				List<String> chrNames = construct.getChromosomesNames();
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
				case BED:case BAM:case GFF:
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

	public static void main(String[] args){
		try {
			ConvertToSQLite c = new ConvertToSQLite("/Users/jarosz/Documents/epfl/flat_files/Rip140_day0_May_treat_afterfiting_chr14.wig",Extension.WIG,70);
			c.convert("/Users/jarosz/Documents/epfl/flat_files/Rip140_day0_May_treat_afterfiting_chr14.sql","qualitative");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ParsingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ExtensionNotRecognisedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
