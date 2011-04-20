package ch.epfl.bbcf.bbcfutils.conversion.json;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.json.JSONArray;
import org.json.JSONException;

import ch.epfl.bbcf.bbcfutils.Utility;
import ch.epfl.bbcf.bbcfutils.conversion.json.pojo.TrackData;
import ch.epfl.bbcf.bbcfutils.conversion.json.pojo.TrackDataQuantitative;
import ch.epfl.bbcf.bbcfutils.conversion.json.pojo.TrackData.ClientConfig;
import ch.epfl.bbcf.bbcfutils.conversion.json.pojo.TrackData.Type;
import ch.epfl.bbcf.bbcfutils.json.JsonMapper;
import ch.epfl.bbcf.bbcfutils.parser.exception.ParsingException;
import ch.epfl.bbcf.bbcfutils.parser.feature.ExtendedQualitativeFeature;
import ch.epfl.bbcf.bbcfutils.parser.feature.JSONFeature;
import ch.epfl.bbcf.bbcfutils.parser.feature.QualitativeFeature;
import ch.epfl.bbcf.bbcfutils.sqlite.SQLiteAccess;


public class ConvertToJSON {

	public enum Extension {GFF,BED,BAM}


	private String inputPath;
	SQLiteAccess access;

	private File curOuput;

	private int curChunkSize;
	private int curChunkNumber;
	private int start;
	private Type processingType;
	private int end;
	private String type;

	public ConvertToJSON(String inputPath,String type){
		this.setInputPath(inputPath);
		this.type = type;
	}


	public boolean convert(String outputPath,String dbName,String ressourceUrl,String trackName) throws ParsingException{
		File test = new File(outputPath+"/"+dbName);
		if(test.exists()){
			System.err.println("already exist");
		}
		boolean success = (new File(outputPath+"/"+dbName)).mkdir();
		if(!success){
			System.err.println("Directory creation failed : "+outputPath+"/"+dbName);
			return false;
		}
		try {
			this.access = (SQLiteAccess) SQLiteAccess.getConnectionWithDatabase(inputPath);
		} catch (InstantiationException e) {
			throw new ParsingException(e);
		} catch (IllegalAccessException e) {
			throw new ParsingException(e);
		} catch (ClassNotFoundException e) {
			throw new ParsingException(e);
		} catch (SQLException e) {
			throw new ParsingException(e);
		}



		//get chromosomes
		Map<String, Integer> chromosomes;
		try {
			chromosomes = access.getChromosomesAndLength();
		} catch (SQLException e1) {
			throw new ParsingException(e1);
		}
		if(null==chromosomes){
			System.err.println("no chromosomes found");
			return false;
		}
		if(type.equalsIgnoreCase("quantitative")){
			try{
				for(Map.Entry<String, Integer> entry : chromosomes.entrySet()){
					final String chromosome = entry.getKey();
					makeDirectory(chromosome,outputPath+"/"+dbName);
					float max;
					max = access.getMaxScoreForChr(chromosome);
					float min = access.getMinScoreForChr(chromosome);
					TrackDataQuantitative tdata = new TrackDataQuantitative(
							Constants.TILE_WIDTH, Math.round(min), Math.round(max), dbName, chromosome);
					String json = JsonMapper.serialize(tdata);
					OutputStream out = new FileOutputStream(curOuput+"/../"+chromosome+".json");
					Utility.write(json,out);
					Utility.close(out);
				}
				access.close();
			} catch (SQLException e) {
				throw new ParsingException(e);
			} catch (IOException e) {
				throw new ParsingException(e);
			}
			return true;
		} else if(type.equalsIgnoreCase("qualitative")){
			//guess type (BASIC or EXTENDED)
			processingType = Type.EXTENDED;
			if(!chromosomes.isEmpty()){
				String firstChromosome = chromosomes.keySet().iterator().next();
				try{
					this.access.getExtendedQualitativeFeature(firstChromosome);
				}catch(SQLException e){
					processingType = Type.BASIC;
				}
			}
			System.out.println("TYPE : "+processingType);


			//iterate throuht chromosomes
			for(Map.Entry<String, Integer> entry : chromosomes.entrySet()){
				final String chromosome = entry.getKey();
				final int length = entry.getValue();
				//make directory
				if(!makeDirectory(chromosome,outputPath+"/"+dbName)){
					System.err.println("cannot create directory : "+outputPath+"/"+dbName+"/"+chromosome);
					return false;
				}
				try {
					JSONChromosome jsonChromosome = new JSONChromosome(chromosome);
					ResultSet r;

					r = access.prepareQualitativeFeatures(chromosome);

					//iterate throught features
					List<String> types = null;
					switch(processingType){
					case BASIC:
						while(r.next()){
							QualitativeFeature feature = access.getNextQualitativeFeature(r);
							JSONFeature jsonFeature = new JSONFeature(feature);
							jsonChromosome.addFeature(jsonFeature);
						}
						break;
					case EXTENDED:
						types = new ArrayList<String>();
						while(r.next()){
							ExtendedQualitativeFeature feature = access.getNextExtendedQualitativeFeature(r);
							JSONFeature jsonFeature = new JSONFeature(feature);
							jsonChromosome.addFeature(jsonFeature);
							if(!types.contains(feature.getType())){
								types.add(feature.getType());
							}
						}
						break;
					}
					r.close();
					writeChromosome(jsonChromosome,length,dbName,ressourceUrl,types,trackName);
				} catch (SQLException e) {
					throw new ParsingException(e);
				} catch (JSONException e) {
					throw new ParsingException(e);
				} catch (InstantiationException e) {
					throw new ParsingException(e);
				} catch (IllegalAccessException e) {
					throw new ParsingException(e);
				} catch (ClassNotFoundException e) {
					throw new ParsingException(e);
				} catch (IOException e) {
					throw new ParsingException(e);
				}
			}
			try {
				access.close();
			} catch (SQLException e) {
				throw new ParsingException(e);
			}
			return true;
		}
		return false;
	}








	private boolean writeChromosome(JSONChromosome jsonChromosome,int length,String dbName,String ressourceUrl,List<String> types,String trackName) throws JSONException, SQLException, InstantiationException, IllegalAccessException, ClassNotFoundException, IOException {
		List<JSONFeat> list = new ArrayList<JSONFeat>();
		for(Map.Entry<String, JSONFeat> entry : jsonChromosome.features.entrySet()){
			JSONFeat f = entry.getValue();
			if(f.name!=null && !f.name.equalsIgnoreCase("")){
				f.parentName = f.name;
			} else {
				f.parentName=entry.getKey();
			}
			list.add(f);
		}
		NCList.sort(list);
		list = NCList.arrange(list);
		TrackData trackData = new TrackData(processingType,types);
		curChunkSize=0;
		curChunkNumber=0;
		start=-1;
		end=-1;
		trackData.setFeatureNCList(new ArrayList<Object>());
		//write lazyOutput & NClist
		LazyOutput lazyOutput = new LazyOutput(curOuput.getAbsolutePath());
		for(int i=0;i<list.size();i++){
			trackData = writeLazyFeatures(list.get(i),lazyOutput,trackData);
		}
		addToTrackData(trackData);
		lazyOutput.writeLazyOutputAndClose();
		trackData.setLazyfeatureUrlTemplate("../"+dbName+"/"+jsonChromosome.chromosome_name+"/lazyfeatures-{chunk}.json");
		//get feature count
		int featureCount;
		featureCount = access.getFeatureCountForChromosome(jsonChromosome.chromosome_name);
		//write histogram meta
		trackData = writeHistogramMeta(trackData,jsonChromosome.chromosome_name,dbName,ressourceUrl,length,featureCount);
		//finalize trackData
		trackData = finalizeTrackData(trackData,trackName,types,featureCount);
		String json = JsonMapper.serialize(trackData);
		OutputStream out = new FileOutputStream(curOuput+"/../"+jsonChromosome.chromosome_name+".json");
		Utility.write(json,out);
		Utility.close(out);
		return true;


	}
	private TrackData finalizeTrackData(TrackData trackData,String trackName,List<String> types,int featureCount) {
		//note : headers already written
		trackData.setKey(trackName);
		trackData.setClassName(Constants.CLASSNAME);
		trackData.setClientConfig(new ClientConfig(processingType, types));
		trackData.setType(Constants.TYPE);
		trackData.setLabel(trackName);
		trackData.setFeatureCount(featureCount);
		return trackData;
	}


	private TrackData writeHistogramMeta(TrackData trackData,String chromosome,String dbName,String ressourceUrl,int length,int featureCount) throws FileNotFoundException, InstantiationException, IllegalAccessException, ClassNotFoundException, SQLException {
		double t = (length*2.5)/featureCount;
		int threshold = 0;
		for (int i : Constants.zooms){
			threshold = i;
			if(i>t){
				break;
			}
		}
		trackData = Histo.writeHistogramMeta(
				trackData, length, threshold, dbName, chromosome, ressourceUrl);

		return Histo.writeHistoFiles(trackData,curOuput, chromosome, length, threshold, inputPath);



	}


	/**
	 * handler.writeValues(
					this.chromosome_name, f.start,f.end,f.parentName,f.strand,f.feature, f.featureCount,finish,0);
	 * @throws FileNotFoundException 
	 */

	private TrackData writeLazyFeatures(JSONFeat feat, LazyOutput lazyOutput, TrackData trackData) throws FileNotFoundException {

		if(start==-1){
			start=feat.start;
		}
		lazyOutput.add(feat.feature);
		curChunkSize+=feat.featureCount;
		if(feat.featureCount==0){
			curChunkSize++;
		}
		end=feat.end;
		if(curChunkSize>=Constants.CHUNK_SIZE){
			//close lazy output
			lazyOutput.writeLazyOutputAndClose();
			//write in trackData
			addToTrackData(trackData);
			//re init values
			curChunkSize=0;
			curChunkNumber++;
			start=-1;
			lazyOutput=new LazyOutput(curOuput.getAbsolutePath());
		}
		return trackData;


	}

	private TrackData addToTrackData(TrackData trackData){
		List<Object> featureNCList = trackData.getFeatureNCList();
		featureNCList.add(TrackData.getChunk(start,end,Integer.toString(curChunkNumber)));
		trackData.setFeatureNCList(featureNCList);
		return trackData;
	}
	private boolean makeDirectory(String chromosome,String outputPath) {
		curOuput = new File(outputPath+"/"+chromosome);
		if(!curOuput.mkdirs()){
			return false;
		}
		return true;
	}


	public void setInputPath(String inputPath) {
		this.inputPath = inputPath;
	}
	public String getInputPath() {
		return inputPath;
	}











	public void setType(Type type) {
		this.processingType = type;
	}


	public Type getType() {
		return processingType;
	}












	private class JSONChromosome{
		private String chromosome_name;
		/**
		 * Map containing all the features in the chromosome
		 * key = the id of the feature, value = a list of feature with the same id
		 */
		private Map<String,JSONFeat> features;
		private int randomId;

		public JSONChromosome(String chromosome){
			features = new HashMap<String,JSONFeat>();
			this.chromosome_name = chromosome;
			randomId = 0;
		}
		/**
		 * add a feature to the chromosome and store them in the Map features.
		 * @param feature - the feature to add
		 * @throws JSONException 
		 */
		public void addFeature(JSONFeature feature) throws JSONException{
			JSONFeat feat = new JSONFeat(
					feature.getName(),feature.getStart(),feature.getEnd(),feature.getStrand(),feature.getScore(),feature.getType());
			String key;
			if(feature.getId()==null || feature.getId().equalsIgnoreCase("")){
				randomId++;
				key = Integer.toString(randomId);
			} else {
				key = feature.getId();
			}
			if(features.containsKey(key)){
				feat.merge(features.get(key));
			} else {
				feat.init();
			}
			features.put(key, feat);
		}
	}






	/**
	 * class wich handler a single feature
	 * @author Yohan Jarosz
	 *
	 */
	public class JSONFeat implements Comparable<JSONFeat>{
		private int start,end,strand,featureCount;
		private String type;
		private String parentName;

		public JSONArray subfeatures;
		public JSONArray feature;
		private String name;
		private float score;


		public JSONFeat(String name, int start2, int end2,int strand,float score,String type) {
			this.start = start2;
			this.end = end2;
			this.strand = strand;
			this.type = type;
			this.featureCount=0;
			subfeatures = new JSONArray();
			this.name = name;
			this.score = score;
		}

		/**
		 * merge the features
		 * for instance, if we have two features with the same id A(start,stop,strand,type) :
		 * A(0,50,1,exon) & A(60,70,1,CDS)
		 * we have to merge and it's become [0,70,A,1,[[0,50,1,exon],[60,70,1,CDS]]]
		 * @param old - the feature to merge with the same id
		 * @throws JSONException
		 */
		public void merge(JSONFeat old) throws JSONException {
			this.subfeatures = old.subfeatures;
			this.featureCount = old.featureCount+1;
			switch(processingType){
			case BASIC:
				this.subfeatures.put(new JSONArray("["+start+","+end+"]"));
				break;
			case EXTENDED :
				this.subfeatures.put(new JSONArray("["+start+","+end+","+strand+",\""+type+"\"]"));
				break;
			}

			if(old.start<this.start){
				this.start=old.start;
			}
			if(old.end>this.end){
				this.end=old.end;
			}
		}
		/**
		 * initialize the feature if there is no existing id for it
		 * @throws JSONException
		 */
		public void init() throws JSONException {
			this.featureCount++;
			switch(processingType){
			case BASIC:
				this.subfeatures.put(new JSONArray("["+start+","+end+"]"));
				break;
			case EXTENDED :
				this.subfeatures.put(new JSONArray("["+start+","+end+","+strand+",\""+type+"\"]"));
				break;
			}
		}





		public void initFeature() throws JSONException{
			switch(processingType){
			case BASIC:
				this.feature=new JSONArray("["+start+","+end+","+score+","+Utility.protect(parentName)+","+strand+"]");
				break;
			case EXTENDED:
				this.feature=new JSONArray("["+start+","+end+","+Utility.protect(parentName)+","+strand+"]");
				if(featureCount>0){
					this.feature.put(this.subfeatures);
				}
				break;
			}


		}


		/**
		 * method to sort the features
		 */
		@Override
		public int compareTo(JSONFeat f) {
			int res = this.start - f.start;
			if(res == 0){
				return f.end - this.end;
			}
			return res;
		}
		public int getStart(){
			return start;
		}
		public int getEnd(){
			return end;
		}
		public int getFeatureCount(){
			return featureCount;
		}
		public void setFeatureCount(int i){
			this.featureCount = i;
		}

		/**
		 * nest a feature into another
		 * @param nextFeat
		 * @throws JSONException
		 */
		public void nesting(JSONFeat nextFeat) throws JSONException {
			JSONArray nested = null;
			if(this.feature.length()==Constants.SUBLIST_INDEX+1){
				nested = this.feature.getJSONArray(Constants.SUBLIST_INDEX-1);
			}
			if(null!=nested){
				nested.put(nextFeat.feature);
			} else {
				nested = new JSONArray();
				nested.put(nextFeat.feature);
			}
			this.feature.put(Constants.SUBLIST_INDEX-1,nested);
			//update count
			this.featureCount+=nextFeat.featureCount;
		}



	}




	private class LazyOutput{
		private final JSONArray array;
		private final String output;

		protected LazyOutput(final String output){
			this.output=output;
			this.array = new JSONArray();
		}
		protected void add(final JSONArray feature){
			this.array.put(feature);
		}
		/**
		 * write the chunk to the output and then
		 * close it
		 * @throws FileNotFoundException 
		 */
		protected void writeLazyOutputAndClose() throws FileNotFoundException{
			final File tmp = new File(output+"/lazyfeatures-"+curChunkNumber+".json");
			final OutputStream out = new FileOutputStream(tmp,true);
			Utility.write(this.array.toString(),out);
			Utility.close(out);
		}
	}

	public static void main(String[] args){
		ConvertToJSON c = new ConvertToJSON("/Users/jarosz/Documents/epfl/flat_files/Rip140_day0_May_treat_afterfiting_chr14.sql","quantitative");
		try {
			c.convert("/Users/jarosz/Documents/epfl/flat_files/gff/json","14JSONZ", "../../tracks_dev", "THE_TRACK");
		} catch (ParsingException e) {
			e.printStackTrace();
		}

	}






}