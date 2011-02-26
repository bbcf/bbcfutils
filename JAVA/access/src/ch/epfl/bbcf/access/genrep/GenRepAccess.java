package ch.epfl.bbcf.access.genrep;

import java.io.IOException;
import java.net.URL;
import java.net.URLConnection;

import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;

import ch.epfl.bbcf.access.utility.SelectOption;



/**
 * class which handle the connections and queries 
 * to Genrep. A database wich stores information aubout genomes
 * for the BBCF
 * @author Yohan Jarosz
 *
 */
public class GenRepAccess implements GenRepKeys{

	protected static final JSONObject no_json_result = null;
	protected static final JSONArray no_json_array_result = null;

	public static boolean testConnection() throws IOException{
		URL url;
		URLConnection urlConnection = null;
		url = new URL(GENEREP+ASSEMBLIES_URL+GET);
		urlConnection = url.openConnection();
		return urlConnection!=null;
	}
	
	protected final static String SHOW(final int id){
		return "/"+id+".json";
	}
	
	public static NR_AssembliesAccess getNR_assemblies() throws JSONException, IOException{
		return new NR_AssembliesAccess();
	}
	public static GenomesAccess getGenomeAccess() throws IOException, JSONException {
		return new GenomesAccess();
	}
	
	
	
	/**
	 * convenient method to get a SelectOption from a JSONObject
	 * @param json - the JSONObject
	 * @param key - the key, it must return an int with the method json.getInt(key) 
	 * @param value - the value, it must return a String with the method json.getString(value)
	 * @return a SelectOption (getInt(key),getString(value))
	 * @throws JSONException
	 */
	public SelectOption getSelectOption(JSONObject json,String key,String value) throws JSONException{
		return new SelectOption(json.getInt(key),json.getString(value));
	}

	
	
	
	
	
	
	
	
}
