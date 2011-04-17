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
import ch.epfl.bbcf.bbcfutils.conversion.json.pojo.TrackData.ClientConfig;
import ch.epfl.bbcf.bbcfutils.conversion.json.pojo.TrackData.Type;
import ch.epfl.bbcf.bbcfutils.json.JsonMapper;
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

	public ConvertToJSON(String inputPath,Type type){
		this.setInputPath(inputPath);
		this.setType(type);
	}


	public boolean convert(String outputPath,String dbName,String ressourceUrl,String trackName){
		File dir = new File(outputPath+"/"+dbName);
		if(!dir.mkdir()){
			System.err.println("Directory creation failed : "+dir.getAbsolutePath());
			return false;
		}
		try {
			this.access = (SQLiteAccess) SQLiteAccess.getConnectionWithDatabase(inputPath);
		} catch (InstantiationException e) {
			e.printStackTrace();
			return false;
		} catch (IllegalAccessException e) {
			e.printStackTrace();
			return false;
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
			return false;
		} catch (SQLException e) {
			e.printStackTrace();
			return false;
		}
		//get chromosomes
		Map<String, Integer> chromosomes;
		try {
			chromosomes = access.getChromosomesAndLength();
		} catch (SQLException e) {
			e.printStackTrace();
			return false;
		}
		if(null==chromosomes){
			System.err.println("no chromosomes found");
			return false;
		}
		//iterate throuht chromosomes
		for(Map.Entry<String, Integer> entry : chromosomes.entrySet()){
			final String chromosome = entry.getKey();
			final int length = entry.getValue();
			//make directory
			if(!makeDirectory(chromosome,outputPath+"/"+dbName)){
				System.err.println("cannot create directory : "+outputPath+"/"+dbName+"/"+chromosome);
				return false;
			}

			JSONChromosome jsonChromosome = new JSONChromosome(chromosome);
			try {
				ResultSet r = access.prepareQualitativeFeatures(chromosome);
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
				e.printStackTrace();
			} catch (JSONException e) {
				e.printStackTrace();
			}




		}

		try {
			access.close();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return true;
	}








	private boolean writeChromosome(JSONChromosome jsonChromosome,int length,String dbName,String ressourceUrl,List<String> types,String trackName) {
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
		try {
			list = NCList.arrange(list);
			
		} catch (JSONException e) {
			e.printStackTrace();
		}

		TrackData trackData = new TrackData(processingType,types);
		curChunkSize=0;
		curChunkNumber=0;
		start=-1;
		trackData.setFeatureNCList(new ArrayList<Object>());
		//write lazyOutput & NClist
		LazyOutput lazyOutput = new LazyOutput(curOuput.getAbsolutePath());
		for(int i=0;i<list.size();i++){
			trackData = writeLazyFeatures(list.get(i),lazyOutput,trackData);
		}
		lazyOutput.writeLazyOutputAndClose();
		trackData.setLazyfeatureUrlTemplate("../"+dbName+"/"+jsonChromosome.chromosome_name+"/lazyfeatures-{chunk}.json");
		//get feature count
		int featureCount;
		try {
			featureCount = access.getFeatureCountForChromosome(jsonChromosome.chromosome_name);
		} catch (SQLException e) {
			e.printStackTrace();
			return false;
		}
		//write histogram meta
		trackData = writeHistogramMeta(trackData,jsonChromosome.chromosome_name,dbName,ressourceUrl,length,featureCount);
		//finalize trackData
		trackData = finalizeTrackData(trackData,trackName,types);
		try {
			String json = JsonMapper.serialize(trackData);
			OutputStream out = new FileOutputStream(curOuput+"/"+jsonChromosome.chromosome_name+".json");
			Utility.write(json,out);
			Utility.close(out);
		} catch (IOException e) {
			e.printStackTrace();
		}
		return true;
		

	}
	private TrackData finalizeTrackData(TrackData trackData,String trackName,List<String> types) {
		//note : headers already written
		trackData.setKey(trackName);
		trackData.setClassName(Constants.CLASSNAME);
		trackData.setClientConfig(new ClientConfig(processingType, types));
		trackData.setType(Constants.TYPE);
		trackData.setLabel(trackName);
		return trackData;
	}


	private TrackData writeHistogramMeta(TrackData trackData,String chromosome,String dbName,String ressourceUrl,int length,int featureCount) {
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
		
		try {
			return Histo.writeHistoFiles(trackData,curOuput, chromosome, length, threshold, inputPath);
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (InstantiationException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		} catch (SQLException e) {
			e.printStackTrace();
		}
		return trackData;
		
		
		
	}


	/**
	 * handler.writeValues(
					this.chromosome_name, f.start,f.end,f.parentName,f.strand,f.feature, f.featureCount,finish,0);
	 */

	private TrackData writeLazyFeatures(JSONFeat feat, LazyOutput lazyOutput, TrackData trackData) {

		List<Object> featureNCList = trackData.getFeatureNCList();
		if(start==-1){
			start=feat.start;
		}
		//TODO
		lazyOutput.add(feat.feature);
		curChunkSize+=feat.featureCount;
		if(feat.featureCount==0){
			curChunkSize++;
		}
		if(curChunkSize>=Constants.CHUNK_SIZE){
			//close lazy output
			lazyOutput.writeLazyOutputAndClose();
			//write in trackData
			featureNCList.add(TrackData.getChunk(start,feat.end,Integer.toString(curChunkNumber)));
			trackData.setFeatureNCList(featureNCList);
			//re init values
			curChunkSize=0;
			curChunkNumber++;
			start=-1;
			lazyOutput=new LazyOutput(curOuput.getAbsolutePath());
		}
		return trackData;


	}


	private boolean makeDirectory(String chromosome,String outputPath) {
		curOuput = new File(outputPath+"/"+chromosome);
		if(!curOuput.mkdir()){
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
			System.out.println("merge");
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
		 */
		protected void writeLazyOutputAndClose(){
			final File tmp = new File(output+"/lazyfeatures-"+curChunkNumber+".json");
			try {
				final OutputStream out = new FileOutputStream(tmp,true);
				Utility.write(this.array.toString(),out);
				Utility.close(out);
			} catch (final FileNotFoundException e) {
				e.printStackTrace();
			}
		}
	}

	public static void main(String[] args){
		ConvertToJSON c = new ConvertToJSON("/Users/jarosz/Documents/epfl/flat_files/gff/Mus_musculus.sql", Type.BASIC);
		c.convert("/Users/jarosz/Documents/epfl/flat_files/gff/json","be757fd10542050d7ed9fd29edcd9200.db", "../../tracks_dev", "THE_TRACK");
	
	}


	



}