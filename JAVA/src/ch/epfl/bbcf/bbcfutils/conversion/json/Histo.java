package ch.epfl.bbcf.bbcfutils.conversion.json;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.OutputStream;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

import ch.epfl.bbcf.bbcfutils.Utility;
import ch.epfl.bbcf.bbcfutils.conversion.json.pojo.TrackData;
import ch.epfl.bbcf.bbcfutils.conversion.json.pojo.TrackData.ArrayParam;
import ch.epfl.bbcf.bbcfutils.conversion.json.pojo.TrackData.HistStats;
import ch.epfl.bbcf.bbcfutils.conversion.json.pojo.TrackData.HistogramMeta;
import ch.epfl.bbcf.bbcfutils.sqlite.SQLiteAccess;


public class Histo {



	public static TrackData writeHistogramMeta(TrackData trackData,int chrLength, int threshold,
			String outputName, String chr,String ressourceUrl) {
		List<TrackData.HistogramMeta> histogramMeta = new ArrayList<TrackData.HistogramMeta>();
		int i = (int) Math.ceil(chrLength/threshold); 
		String url = ressourceUrl+"/"+outputName+"/"+chr+"/hist-"+threshold+"-{chunk}.json";
		ArrayParam ar = new ArrayParam(i,10000,url);
		HistogramMeta hm = new HistogramMeta(threshold, ar);
		histogramMeta.add(hm);
		//		int megathreshold = threshold*100;
		//		int im = (int) Math.ceil(chrLength/megathreshold); 
		//		String urlm = ressourceUrl+"/"+outputName+"/"+chr+"/hist-"+megathreshold+"-{chunk}.json";
		//		ArrayParam arm = new ArrayParam(im,10000,urlm);
		//		HistogramMeta hmm = new HistogramMeta(threshold, arm);
		//		histogramMeta.add(hmm);
		trackData.setHistogramMeta(histogramMeta);
		return trackData;
	}







	public static TrackData writeHistoFiles(TrackData trackData,File outDir,String chromosome,
			int chrLength,int threshold,String databasePath) throws FileNotFoundException, InstantiationException, IllegalAccessException, ClassNotFoundException, SQLException {


		List<Integer> bigArray = new ArrayList<Integer>();
		//init ouputs
		int curChunk = 0;
		int curPos = 0;
		OutputStream outT = newHistChunkOutPut(outDir,threshold,curChunk);
		List<Integer> counts = new ArrayList<Integer>();
		int megathreshold = threshold*100;
		int mCurPos = 0;
		int mCurChunk = 0;
		int mCount = 0;
		OutputStream outMT = newHistChunkOutPut(outDir,megathreshold,mCurChunk);
		List<Integer> mCounts = new ArrayList<Integer>();
		//init connections
		SQLiteAccess access = (SQLiteAccess) SQLiteAccess.getConnectionWithDatabase(databasePath);
		ResultSet r = access.getStartEndForChromosome(chromosome);
		List<Integer> starts = new ArrayList<Integer>();
		List<Integer> ends = new ArrayList<Integer>();
		while(r.next()){
			starts.add(r.getInt(1));
			ends.add(r.getInt(2));
		}
		r.close();
		access.close();
		int startCursor = 0;
		int endCursor = 0;

		for(int i=1;i<=chrLength;i+=threshold){
			int start = i;
			int end = i+threshold;
			startCursor = moveStartcursor(start,end,startCursor,starts,ends);
			endCursor = moveEndCursor(start,end,startCursor,endCursor,starts,ends);
			int count = getCount(start,end,startCursor,endCursor,starts,ends);
			counts.add(count);
			bigArray.add(count);
			curPos++;
			mCount+=count;
			if(curPos%100==0){
				mCounts.add(mCount);
				mCount=0;
				mCurPos++;
			}
			if(curPos%Constants.HIST_CHUNK_SIZE==0){
				Utility.write(counts.toString(),outT);
				counts=new ArrayList<Integer>();
				Utility.close(outT);
				curChunk++;
				outT = newHistChunkOutPut(outDir,threshold,curChunk);
			}
			if(mCurPos%Constants.HIST_CHUNK_SIZE==0 && mCurPos!=0){
				Utility.write(mCounts.toString(),outMT);
				mCounts=new ArrayList<Integer>();
				Utility.close(outMT);
				mCurChunk++;
				outMT = newHistChunkOutPut(outDir,megathreshold,mCurChunk);
			}
		}
		//finalize
		if(outT!=null){
			Utility.write(counts.toString(),outT);
			Utility.close(outT);
		}
		if(outMT!=null){
			Utility.write(mCounts.toString(),outMT);
			Utility.close(outMT);
		}
		return writeHistoStats(trackData,threshold,bigArray,chrLength);


	}

	protected static int getCount(int start, int end, int startCursor,
			int endCursor, List<Integer> starts, List<Integer> ends) {
		if((starts.get(startCursor)<=end && ends.get(startCursor)>=start)){
			return endCursor -startCursor +1;
		}
		return 0;
	}



	protected static int moveEndCursor(int start, int end, int startCursor,int endCursor, List<Integer> starts,List<Integer> ends) {
		for(int i=endCursor;i<starts.size();i++){
			endCursor=i;
			if(starts.get(i)>end){
				if(endCursor==0){
					return endCursor;
				}
				return endCursor-1;
			}
		}
		return endCursor;	
	}



	protected static int moveStartcursor(int start, int end, int startCursor, List<Integer> starts,List<Integer> ends) {
		for(int i=startCursor;i<ends.size();i++){
			startCursor=i;
			if(ends.get(i)>=start){
				return startCursor;
			}
		}
		return startCursor;
	}




	protected static TrackData writeHistoStats(TrackData trackData,int threshold, List<Integer> bigArray, int chrLength) {
		List<HistStats> stats = new ArrayList<HistStats>();
		for(int z :Constants.zooms){
			int base = z*threshold;
			if(base<chrLength){
				HistStats stat = new HistStats();
				List<Integer>baseArray= new ArrayList<Integer>();
				for(int i=0;i<bigArray.size()-z-1;i+=z){
					List<Integer>tmpArray=bigArray.subList(i, i+z);
					int c=0;
					for(int j :tmpArray){
						c+=j;
					}
					baseArray.add(c);
				}
				int max = 0;
				int mean = 0;
				for(int k:baseArray){
					mean+=k;
					if(k>max){
						max = k;
					}
				}
				try{
					mean = mean/baseArray.size();
					stat.setBases(base);
					stat.setMax(max);
					stat.setMean(mean);
					stats.add(stat);
				} catch(ArithmeticException e){}
			} else {
				break;
			}
		}
		trackData.setHistStats(stats);
		return trackData;
	}



	protected static OutputStream newHistChunkOutPut(File outputDir,
			int threshold, int curChunk) throws FileNotFoundException {

		return new FileOutputStream(
				new File(outputDir.getAbsolutePath()+"/hist-"+threshold+"-"+curChunk+".json"));
	}





	/**
	 * write the json output
	 * @param name2 - name of the track
	 * @param chrLength - length of the chromosome
	 * @param threshold - threshold
	 * @param featureCount - nb of features
	 * @param inputFilePath - the input database path
	 * @param fileName - the file name
	 * @param outputPath - the ouput directory
	 * @param ressourceUrl -the ressource url for the browser
	 * @param chr - the chromosome name
	 * @param out - an ouputsttream
	 * @param clientConfig2 - the clientconfig
	 * @throws FileNotFoundException
	 * @throws InstantiationException
	 * @throws IllegalAccessException
	 * @throws ClassNotFoundException
	 * @throws SQLException
	 */
}
