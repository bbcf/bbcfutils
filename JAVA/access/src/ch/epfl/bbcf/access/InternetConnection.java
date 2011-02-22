package ch.epfl.bbcf.access;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLConnection;


/**
 * class which handle POST and GET connections
 * @author Yohan Jarosz
 *
 */
public class InternetConnection {

	public final static String MIME_TYPE_TEXT = "text/plain";
	public final static String MIME_TYPE_FORM_APPLICATION = "application/x-www-form-urlencoded";
	/**
	 * Test an URL connection
	 * @param adress - the URI
	 * @return
	 * @throws IOException
	 */
	public static boolean testConnection(String adress) throws IOException{
		URL url;
		URLConnection urlConnection = null;
		url = new URL(adress);
		urlConnection = url.openConnection();
		return urlConnection!=null;
	}

	/**
	 * Open an URL connection with the method provided
	 * @param adress - the URI
	 * @param method - the method ("POST","GET","PUSH",...)
	 * @return an HTTP connection
	 * @throws IOException
	 */
	private static HttpURLConnection openConnection(String adress,String method) throws IOException{
		URL url = new URL(adress);
		HttpURLConnection urlConnection = (HttpURLConnection) url.openConnection();
		urlConnection.setRequestMethod(method);
		return urlConnection;
	}

	/**
	 * get back the stream of the connection : can be the normal input or an error stream
	 * @param urlConnection
	 * @return
	 */
	private static DataInputStream getConnectionStream(HttpURLConnection urlConnection){
		DataInputStream inStream = null;
		try {
			inStream  = new DataInputStream(urlConnection.getInputStream());
		} catch (IOException e) {
			inStream = new DataInputStream(urlConnection.getErrorStream());
		} finally {
			return inStream;
		}

	}

	/**
	 * Convert a stream to a String
	 * @param stream
	 * @return
	 * @throws IOException 
	 */
	private static String convertStreamToString(DataInputStream stream) throws IOException{
		String buffer,result = "";
		InputStreamReader isr = new InputStreamReader(stream);
		BufferedReader br = new BufferedReader(isr);
		while(null!=(buffer = br.readLine())){
			result+=buffer+"\n";
		}
		stream.close();
		return result;
	}
	
	public static String sendGETConnection(final String adress) throws IOException{
		if(testConnection(adress)){
			HttpURLConnection urlConnection = openConnection(adress,"GET");
			DataInputStream inStream = getConnectionStream(urlConnection);
			String result = convertStreamToString(inStream);
			return result;
		}
		return "no connection to "+adress;
	}
	/**
	 * write a body in a URLConnection
	 * @param urlConnection - the connection
	 * @param body - the body
	 * @param mimeType - the Mimetype of the body
	 * @throws IOException
	 */
	private static void writeBody(HttpURLConnection urlConnection,String body,String mimeType) throws IOException{
		System.out.println(body);
		urlConnection.setDoInput(true);
		urlConnection.setDoOutput(true);
		urlConnection.setUseCaches(false);
		urlConnection.setRequestProperty("Content-Type",mimeType);
		urlConnection.setRequestProperty("Content-Length",Integer.toString(body.length()));
		DataOutputStream outStream = new DataOutputStream(urlConnection.getOutputStream());
		outStream.writeBytes(body);
		outStream.flush();
		outStream.close();
	}
	
	/**
	 * Send a POST request
	 * @param adress - the URI
	 * @param body - the body to write
	 * @param mimeType - the Mimetypee of the body
	 * @return
	 * @throws IOException
	 */
	public static String sendPOSTConnection(final String adress,final String body,final String mimeType) throws IOException{
		if(testConnection(adress)){
			HttpURLConnection urlConnection = openConnection(adress,"POST");
			writeBody(urlConnection, body, mimeType);
			DataInputStream inStream = getConnectionStream(urlConnection);
			String result = convertStreamToString(inStream);
			return result;
		}
		return "no connection to "+adress;
	}

	public static void main(String[] args){
		System.out.println("TEST");
		String body = "id=78&data=myData";
		try {
			System.out.println(InternetConnection.sendPOSTConnection("http://myServer.com/post",body,MIME_TYPE_FORM_APPLICATION));
		} catch (IOException e) {
			e.printStackTrace();
		}

	}
}
