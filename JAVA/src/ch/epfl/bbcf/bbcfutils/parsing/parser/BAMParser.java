package ch.epfl.bbcf.bbcfutils.parsing.parser;
//package ch.epfl.bbcf.bbcfutils.parser;
//
//
//import java.io.File;
//import java.io.IOException;
//import java.util.Iterator;
//
//
//import ch.epfl.bbcf.bbcfutils.parser.Handler;
//import ch.epfl.bbcf.bbcfutils.parser.Parser;
//import ch.epfl.bbcf.bbcfutils.parser.exception.ParsingException;
//import ch.epfl.bbcf.bbcfutils.parser.feature.BAMFeature;
//import ch.epfl.bbcf.bbcfutils.parser.feature.Track;
//
//import net.sf.samtools.SAMFileHeader;
//import net.sf.samtools.SAMFileReader;
//import net.sf.samtools.SAMFileSpan;
//import net.sf.samtools.SAMRecord;
//import net.sf.samtools.SAMRecordIterator;
//import net.sf.samtools.util.RuntimeIOException;
//
//
//public class BAMParser extends Parser{
//
//	private Track cur_track;
//
//	public BAMParser(Processing type) {
//		super(type);
//	}
//	@Override
//	protected void processLine(String line, Handler handler)
//	throws ParsingException {
//		//not used
//	}
//	@Override
//	public void parse(File inputSamOrBamFile,Handler handler) throws IOException, ParsingException,RuntimeIOException{
////		final SAMFileReader inputSam_ = new SAMFileReader(inputSamOrBamFile);
////		int cpt = 0;
////		if(inputSam_!=null){
////			SAMRecordIterator it_ = inputSam_.iterator();
////			while(it_.hasNext()){
////				it_.next();
////				cpt++;
////			}
////			inputSam_.close();
////		}
////		int cur = 0;
////		int percent = 0;
////		int coef = cpt/20;
////		
////		final SAMFileReader inputSam = new SAMFileReader(inputSamOrBamFile);
////		if(inputSam!=null){
////			handler.start();
////			cur_track = new Track();
////			newTrack(handler, cur_track);
////			SAMRecordIterator it = inputSam.iterator();
////			while(it.hasNext()){
////				SAMRecord samRecord = it.next();
////				cur++;
////				if(cur%coef==0){
////					percent+=5;
////				}
////				BAMFeature feature = new BAMFeature();
////				feature.setReadName(samRecord.getReadName());
////				feature.setRefName(samRecord.getReferenceName());
////				feature.setStart(samRecord.getAlignmentStart());
////				feature.setStop(samRecord.getAlignmentEnd());
////				newFeature(handler, cur_track, feature);
////			}
////			inputSam.close();
////			handler.end();
////		}
//	}
//
//	public static void main(String[]args){
////		System.out.println("START");
////		final SAMFileReader inputSam = new SAMFileReader(
////				new File("/Users/jarosz/Documents/epfl/flat_files/Ste12_chrII_fwd.bam"));
////		System.out.println("file opened");
////		SAMRecordIterator it = inputSam.iterator();
////		int cpt=0;
////		while(it.hasNext()){
////			final SAMRecord rec = it.next();
////			cpt++;
////		}
////		inputSam.close();
////		
////		
////		SAMFileReader inputSam2 = new SAMFileReader(
////				new File("/Users/jarosz/Documents/epfl/flat_files/Ste12_chrII_fwd.bam"));
////		System.out.println("again");
////		SAMRecordIterator it2 = inputSam2.iterator();
////		int cur = 0;
////		int percent = 0;
////		int coef = cpt/20;
////		while(it2.hasNext()){
////			SAMRecord rec = it2.next();
////			cur++;
////			if(cur%coef==0){
////				percent+=5;
////				System.out.println(percent+"%");
////			}
//////			
////		}
////		inputSam2.close();
//	}
//}
