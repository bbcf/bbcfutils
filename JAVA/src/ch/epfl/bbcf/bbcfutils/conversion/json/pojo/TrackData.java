package ch.epfl.bbcf.bbcfutils.conversion.json.pojo;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;


import ch.epfl.bbcf.bbcfutils.conversion.json.Constants;


public class TrackData {

	private List<String> headers;

	private List<String> subfeatureHeaders;
	public final int subListIndex = Constants.SUBLIST_INDEX;
	private final int lazyIndex = 2;

	private List<HistogramMeta> histogramMeta;
	private List<HistStats> histStats;
	private Map<String,String> subfeatureClasses;
	private List<Object> featureNCList;
	private ClientConfig clientConfig;
	private int featureCount;
	private String key,arrowheadClass,type,
	label,lazyfeatureUrlTemplate,className;

	public enum Type {BASIC,EXTENDED}


	public TrackData(Type type,List<String> types){
		List<String>heads = new ArrayList<String>();
		switch(type){
		case BASIC:
			heads.add("start");
			heads.add("end");
			heads.add("score");
			heads.add("name");
			heads.add("strand");

			break;
		case EXTENDED:
			heads.add("start");
			heads.add("end");
			heads.add("name");
			heads.add("strand");
			heads.add("subfeatures");
			List<String>sheads = new ArrayList<String>();
			sheads.add("start");
			sheads.add("end");
			sheads.add("strand");
			sheads.add("type");
			this.setSubfeatureHeaders(sheads);
			break;
		}
		this.setHeaders(heads);
		
		this.arrowheadClass="transcript-arrowhead";


//		List<TrackData.HistogramMeta> histm = new ArrayList<TrackData.HistogramMeta>();
//		ArrayParam arparam = new ArrayParam(500000, 1000, "../../tracks_dev/be757fd10542050d7ed9fd29edcd9200.db/" +
//		"chr2/hist-500000-{chunk}.json");
//		HistogramMeta hm = new HistogramMeta(500000, arparam);
//		histm.add(hm);
//		this.histogramMeta = histm;
//
//
//
//
//
//		this.histStats = new ArrayList<TrackData.HistStats>();
//		Map<String,String> mabsub = new HashMap<String, String>();
//		mabsub.put("gene", "gene");
//		mabsub.put("cds", "cds");
//		this.setSubfeatureClasses(mabsub);
//		featureNCList=new ArrayList<Object>();
//		featureNCList.add(getChunk(101,200,"3"));
//		featureNCList.add(getChunk(201,300,"4"));
//
//
//
//		this.clientConfig=new ClientConfig();
//		this.featureCount=3000;
//		this.key="testKey";
//		this.arrowheadClass="testarrowheadClass";
//		this.type="testtype";
//		this.label="testlabel";
//		this.lazyfeatureUrlTemplate="testlazyfeatureUrlTemplate";
//		this.className="testclassName";
	}



	/**
	 * Convenient method to build the list of Chunks
	 * @param start - the start 
	 * @param end - the end
	 * @param chunkNumber - the number of the chunk
	 * @return
	 */
	public static List<Object> getChunk(int start,int end,String chunkNumber){
		List<Object> chunk = new ArrayList<Object>();
		chunk.add(start);
		chunk.add(end);
		chunk.add(new Chunk(chunkNumber));
		return chunk;
	}


	/**
	 * The client configuration :
	 * it will contains all about css
	 * @author Yohan Jarosz
	 *
	 */
	public static class ClientConfig{
		private final int labelScale = 5;
		private final int subfeatureScale=10;
		private String featureCallback;
		private String featureCss;

		public ClientConfig(Type type,List<String> types) {
			switch(type){
			case BASIC :
				featureCss = 
					"\"background-color: #0000FF; " +
					"height: 5px;\",\"histCss\":\"background-color: #0000FF;\"";
			case EXTENDED :
				featureCallback =
					"(function(feat, fields, div) { " +
					"if (fields.type){" +
					"getFeatureStyle(feat[fields.type],div);" +
					"}}" +
					")";
			}
		}
		public int getLabelScale() {
			return labelScale;
		}
		public int getSubfeatureScale() {
			return subfeatureScale;
		}
		public void setFeatureCallback(String featureCallback) {
			this.featureCallback = featureCallback;
		}
		public String getFeatureCallback() {
			return featureCallback;
		}
		public void setFeatureCss(String featureCss) {
			this.featureCss = featureCss;
		}
		public String getFeatureCss() {
			return featureCss;
		}

	}

	/**
	 * Some parameters about the histograms
	 * @author Yohan Jarosz
	 *
	 */
	public static class ArrayParam{
		private int length,chunkSize;
		private String urlTemplate;

		public ArrayParam(int l,int cs,String url){
			this.length=l;
			this.chunkSize=cs;
			this.urlTemplate=url;
		}
		public void setLength(int length) {
			this.length = length;
		}
		public int getLength() {
			return length;
		}
		public void setChunkSize(int chunkSize) {
			this.chunkSize = chunkSize;
		}
		public int getChunkSize() {
			return chunkSize;
		}
		public void setUrlTemplate(String urlTemplate) {
			this.urlTemplate = urlTemplate;
		}
		public String getUrlTemplate() {
			return urlTemplate;
		}
	}
	/**
	 * Some parameters about the histograms
	 * @author Yohan Jarosz
	 *
	 */
	public static class HistogramMeta {

		public HistogramMeta(int bpb,ArrayParam a){
			this.basesPerBin=bpb;
			this.arrayParams=a;
		}

		private int basesPerBin;
		private ArrayParam arrayParams;

		public void setBasesPerBin(int basesPerBin) {
			this.basesPerBin = basesPerBin;
		}
		public int getBasesPerBin() {
			return basesPerBin;
		}
		public void setArrayParams(ArrayParam arrayParams) {
			this.arrayParams = arrayParams;
		}
		public ArrayParam getArrayParams() {
			return arrayParams;
		}


	}

	/**
	 * Some statistics about histograms
	 * @author Yohan Jarosz
	 *
	 */
	public static class HistStats {
		private int bases,max,mean;

		public void setMean(int mean) {
			this.mean = mean;
		}

		public int getMean() {
			return mean;
		}

		public void setBases(int bases) {
			this.bases = bases;
		}

		public int getBases() {
			return bases;
		}

		public void setMax(int max) {
			this.max = max;
		}

		public int getMax() {
			return max;
		}
	}


	/**
	 * a chunk
	 * @author Yohan Jarosz
	 *
	 */
	public static class Chunk {
		public Chunk(String number){
			this.chunk=number;
		}
		private String chunk;

		public void setChunk(String chunk) {
			this.chunk = chunk;
		}

		public String getChunk() {
			return chunk;
		}
	}



	


	//Getter and setters


	public void setHistogramMeta(List<HistogramMeta> histogramMeta) {
		this.histogramMeta = histogramMeta;
	}
	public List<HistogramMeta> getHistogramMeta() {
		return histogramMeta;
	}
	public void setHistStats(List<HistStats> histStats) {
		this.histStats = histStats;
	}
	public List<HistStats> getHistStats() {
		return histStats;
	}
	public void setFeatureCount(int featureCount) {
		this.featureCount = featureCount;
	}
	public int getFeatureCount() {
		return featureCount;
	}
	public void setArrowheadClass(String arrowheadClass) {
		this.arrowheadClass = arrowheadClass;
	}
	public String getArrowheadClass() {
		return arrowheadClass;
	}
	public void setType(String type) {
		this.type = type;
	}
	public String getType() {
		return type;
	}
	public void setKey(String key) {
		this.key = key;
	}
	public String getKey() {
		return key;
	}
	public void setClassName(String className) {
		this.className = className;
	}
	public String getClassName() {
		return className;
	}
	public void setLazyfeatureUrlTemplate(String lazyfeatureUrlTemplate) {
		this.lazyfeatureUrlTemplate = lazyfeatureUrlTemplate;
	}
	public String getLazyfeatureUrlTemplate() {
		return lazyfeatureUrlTemplate;
	}
	public void setLabel(String label) {
		this.label = label;
	}
	public String getLabel() {
		return label;
	}
	public void setClientConfig(ClientConfig clientConfig) {
		this.clientConfig = clientConfig;
	}
	public ClientConfig getClientConfig() {
		return clientConfig;
	}
	public void setHeaders(List<String> headers) {
		this.headers = headers;
	}
	public List<String> getHeaders() {
		return headers;
	}
	public void setSubfeatureHeaders(List<String> subfeatureHeaders) {
		this.subfeatureHeaders = subfeatureHeaders;
	}
	public List<String> getSubfeatureHeaders() {
		return subfeatureHeaders;
	}
	public void setSubfeatureClasses(Map<String,String> subfeatureClasses) {
		this.subfeatureClasses = subfeatureClasses;
	}
	public Map<String,String> getSubfeatureClasses() {
		return subfeatureClasses;
	}
	public int getLazyIndex() {
		return lazyIndex;
	}
	public int getSubListIndex() {
		return subListIndex;
	}
	public void setFeatureNCList(List<Object> featureNCList) {
		this.featureNCList = featureNCList;
	}
	public List<Object> getFeatureNCList() {
		return featureNCList;
	}
}
