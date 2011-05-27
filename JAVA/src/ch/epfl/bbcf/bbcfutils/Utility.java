package ch.epfl.bbcf.bbcfutils;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.math.BigInteger;
import java.security.DigestInputStream;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;

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

	public static void read(InputStream in){
		while (true) {
			synchronized (buffer) {
				int amountRead = 0;
				try {
					amountRead = in.read(buffer);
				} catch (IOException e) {
					e.printStackTrace();
				}
				if (amountRead == -1) {
					break;
				}
			}
		}
		try {
			in.close();
		} catch (IOException e) {
			e.printStackTrace();
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


	public static String protect(String str){
		if(!str.startsWith("\"")){
			return "\""+str+"\"";
		} else {
			return str;
		}
	}


	/**
	 * Get the file digest from a file
	 * @param filePath - the path to the file
	 * @param crypto - the crypto to use ('MD5', 'SHA')
	 * @return the signature in hexa
	 * @throws NoSuchAlgorithmException
	 * @throws IOException
	 */
	public static String getFileDigest(String filePath,String crypto) throws NoSuchAlgorithmException, IOException{
		MessageDigest md = null;
		md = MessageDigest.getInstance(crypto);
		InputStream is = new FileInputStream(filePath);
		try {
			is = new DigestInputStream(is, md);
			read(is);
		}
		finally {
			is.close();
		}
		byte[] digest = md.digest();
		String signature = new BigInteger(1,digest).toString(16);
		return signature;
	}




	public static String getFileMd5(String filePath) throws IOException{
		MessageDigest md = null;
		try {
			md = MessageDigest.getInstance("MD5");
		} catch (NoSuchAlgorithmException e) {
			e.printStackTrace();
		}
		InputStream is = new FileInputStream(filePath);
		try {
			is = new DigestInputStream(is, md);
			read(is);
		}
		finally {
			is.close();
		}
		byte[] digest = md.digest();
		String signature = new BigInteger(1,digest).toString(16);
		return signature;
	}
	public static String getFileMd5(File file) throws IOException{
		return getFileMd5(file.getAbsolutePath());
	}

	public static void main(String[] args){
		try {
			System.out.println(getFileMd5("/Users/jarosz/Desktop/ajax-loader.gif"));
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


}
