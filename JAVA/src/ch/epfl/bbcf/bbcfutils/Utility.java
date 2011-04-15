package ch.epfl.bbcf.bbcfutils;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

public final class Utility {

	protected static final int BUFF_SIZE = 1024;
	protected static final byte[] buffer = new byte[BUFF_SIZE];
	
	/**
	 * convenient method for writing a string to a file
	 * @param toWrite
	 * @param out
	 */
	public static void write(String toWrite,OutputStream out){
		try{
			InputStream in = new ByteArrayInputStream(toWrite.getBytes("UTF-8"));
			while (true) {
				synchronized (buffer) {
					int amountRead = in.read(buffer);
					if (amountRead == -1) {
						break;
					}
					out.write(buffer, 0, amountRead); 
				}
			} 
		}catch (IOException e){
			e.printStackTrace();
		}
	}

	/**
	 * convenient method for copying a strem to another
	 * @param in
	 * @param out
	 */
	protected static void copy(InputStream in, OutputStream out) {
		try {
			while (true) {
				synchronized (buffer) {
					int amountRead = in.read(buffer);
					if (amountRead == -1) {
						break;
					}
					out.write(buffer, 0, amountRead); 
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}finally {
			if (in != null) {
				try {
					in.close();
				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}

	}

	/**
	 * close an outputstream
	 * @param out
	 */
	public static void close(OutputStream out){
		if(out!=null){
			try {
				out.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}

}
