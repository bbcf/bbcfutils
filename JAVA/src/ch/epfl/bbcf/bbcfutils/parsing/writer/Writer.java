package ch.epfl.bbcf.bbcfutils.parsing.writer;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import ch.epfl.bbcf.bbcfutils.parsing.feature.Feature;

public class Writer {

	protected static final byte[]lineSeparator = System.getProperty("line.separator").getBytes();
	
	protected OutputStream out; 

	private static final int BUFF_SIZE = 10*2048;
	private static final byte[] buffer = new byte[BUFF_SIZE];
	private static final int TAB_SIZE = 100;

	private byte[][] features;
	private int cpt;

	protected Writer(){
		features=new byte[TAB_SIZE][];
		cpt=0;
	}


	protected void writeFeature(Feature feature) throws IOException{
		/* add the feature in a bytes tab */
		cpt++;
		features[cpt]=feature.display().getBytes("UTF-8");
		/* if cpt = TAB_SIZE write features & reinit bytes */
		if(cpt>=TAB_SIZE){
			cpt=0;
			writeTab();
			features=new byte[TAB_SIZE][];
		}
	}

	/**
	 * Method called at the end of the writing
	 * @throws IOException 
	 */
	protected void end() throws IOException{
		writeTab();
		out.close();
	}

	/**
	 * write the arrays of bytes to 
	 * the outputstream & add a line separator
	 * @throws IOException
	 */
	private void writeTab() throws IOException{
		for(byte[] tab : features){
			if(null!=tab){
				InputStream in = new ByteArrayInputStream(tab);
				while (true) {
					synchronized (buffer) {
						int amountRead = in.read(buffer);
						if (amountRead == -1) {
							break;
						}
						out.write(buffer, 0, amountRead); 
					}
				} 
				out.write(lineSeparator);
				in.close();
			}
			
		}
	}


}
