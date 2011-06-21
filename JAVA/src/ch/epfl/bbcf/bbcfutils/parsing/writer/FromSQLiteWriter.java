package ch.epfl.bbcf.bbcfutils.parsing.writer;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.Map;

import ch.epfl.bbcf.bbcfutils.exception.ExtensionNotRecognisedException;
import ch.epfl.bbcf.bbcfutils.parsing.Extension;
import ch.epfl.bbcf.bbcfutils.parsing.SQLiteExtension;
import ch.epfl.bbcf.bbcfutils.parsing.feature.BEDFeature;
import ch.epfl.bbcf.bbcfutils.parsing.feature.BioSQLiteQualitative;
import ch.epfl.bbcf.bbcfutils.parsing.feature.BioSQLiteQualitativeExt;
import ch.epfl.bbcf.bbcfutils.parsing.feature.BioSQLiteQuantitative;
import ch.epfl.bbcf.bbcfutils.parsing.feature.GFFFeature;
import ch.epfl.bbcf.bbcfutils.parsing.feature.WIGFeature;
import ch.epfl.bbcf.bbcfutils.sqlite.SQLiteAccess;

public class FromSQLiteWriter extends Writer{



	private SQLiteAccess access;

	public FromSQLiteWriter(){
		super();
	}


	public boolean convert(String dbPath,String ouputPath,Extension extension) throws InstantiationException, IllegalAccessException, ClassNotFoundException, SQLException, IOException, ExtensionNotRecognisedException{
		access = SQLiteAccess.getConnectionWithDatabase(dbPath);
		if(null==access){
			return false;
		}
		out = new FileOutputStream(ouputPath);
		if(null==out){
			return false;
		}
		/* get datatype of SQLite file */
		SQLiteExtension ext = access.getDatatype();
		
		/* guess the best extension */
		if(null==extension){
			extension = guessBetterOutput(ext);
		}
		/* write header line */
		writeHeader(extension);

		/* write features */
		Map<String,Integer> map = access.getChromosomesAndLength();
		switch(extension){
		case BAM:case WIG: System.err.println("extension wiggle & BAM not supported for writing");break;
		case BED:
			for(String chr : map.keySet()){
				ResultSet r = access.prepareFeatures(chr);
				while(r.next()){
					BioSQLiteQualitative f = access.getNextQualitativeFeature(r,chr);
					BEDFeature bed = f.toBEDFeature();
					writeFeature(bed);
				}
				r.close();
			}
			break;
		case BEDGRAPH:
			for(String chr : map.keySet()){
				ResultSet r = access.prepareFeatures(chr);
				while(r.next()){
					BioSQLiteQuantitative f = access.getNextQuantitativeFeature(r,chr);
					WIGFeature wig = f.toWIGFeature();
					writeFeature(wig);
				}
				r.close();
			}
			break;
		case GFF:
			for(String chr : map.keySet()){
				ResultSet r = access.prepareFeatures(chr);
				while(r.next()){
					BioSQLiteQualitativeExt f = access.getNextExtendedQualitativeFeature(r,chr);
					GFFFeature gff = f.toGFFFeature();
					writeFeature(gff);
				}
				r.close();
			}
			
			break;
		}
		end();
		return true;

	}







	private void writeHeader(Extension extension) throws SQLException, IOException {
		Map<String,String> atts = access.getAttributes();
		String h = "";
		switch(extension){
		case BEDGRAPH : h+="track type=bedGraph";break;
		case BED : h+="track ";break;
		case WIG : h+="track type=wiggle_0";break;
		case GFF : h+="track";break;
		case BAM : 
			System.err.println("conversion to BAM not supportted");
			break;
		}
		out.write((h+" ").getBytes());
		if(null!=atts){
			for(Map.Entry<String, String>entry : atts.entrySet()){
				if(!entry.getKey().equalsIgnoreCase("type")){
					out.write((entry.getKey()+"="+entry.getValue()+" ").getBytes());
				}
			}
		}
		out.write(lineSeparator);
	}


	/**
	 * guess which output should be 
	 * for this type of file
	 * @return
	 * @throws SQLException 
	 * @throws ExtensionNotRecognisedException 
	 */
	private Extension guessBetterOutput(SQLiteExtension ext) throws SQLException, ExtensionNotRecognisedException {
		switch(ext){
		case QUALITATIVE:return Extension.BED;
		case QUANTITATIVE:return Extension.BEDGRAPH;
		case QUALITATIVE_EXTENDED:return Extension.GFF;
		}
		return null;
	}


	public static void main(String[] args){
		FromSQLiteWriter writer = new FromSQLiteWriter();
		try {
			writer.convert("/Users/jarosz/Desktop/final1.db", "/Users/jarosz/Desktop/final1.bed", null);
		} catch (InstantiationException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (SQLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ExtensionNotRecognisedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
