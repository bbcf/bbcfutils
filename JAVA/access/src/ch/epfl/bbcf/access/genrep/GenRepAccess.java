package ch.epfl.bbcf.access.genrep;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.codehaus.jackson.type.TypeReference;

import com.sun.tools.example.debug.bdi.MethodNotFoundException;

import ch.epfl.bbcf.access.InternetConnection;
import ch.epfl.bbcf.access.genrep.Constants.FORMAT;
import ch.epfl.bbcf.access.genrep.Constants.KEY;
import ch.epfl.bbcf.access.genrep.Constants.METHOD;
import ch.epfl.bbcf.access.genrep.json.JsonMapper;
import ch.epfl.bbcf.access.genrep.pojo.Assembly;
import ch.epfl.bbcf.access.genrep.pojo.Chromosome;
import ch.epfl.bbcf.access.genrep.pojo.GenrepObject;
import ch.epfl.bbcf.access.genrep.pojo.Organism;




/**
 * class which handle the connections and queries 
 * to Genrep. A database which stores information about genomes
 * for the BBCF
 * @author Yohan Jarosz
 *
 */
public class GenRepAccess {


	
	/**
	 * prepare the URL to access to Genrep
	 * depending of the parameters
	 * @return the URL
	 * @throws MethodNotFoundException 
	 */
	private static final String prepareUrl(String servUrl,METHOD method,FORMAT format,KEY key,String query) throws MethodNotFoundException{
		switch(method){
		case INDEX : servUrl +="/"+key+"."+format+"?"+query; 
		break;
		case SHOW : servUrl +="/"+key+"/"+query+"."+format;
		break;
		case ALL : servUrl +="/"+key+"."+format;
		break;
		default :
			throw new MethodNotFoundException(" method must be part of "+METHOD.values());
		}
		return servUrl;
	
	}
	/**
	 * Launch the query to Genrep and get
	 * back an GenrepObject
	 * @param servUrl - the Genrep server URL
	 * @param method - the methods supported by Genrep
	 * @param format - the format supported by Genrep
	 * @param clazz - the class you want to get back (map the deserialization of the JSON)
	 * @param key - the key to put in the URL (usually the same than the className)
	 * @param query - the query (e.g : an id,md5,...)
	 * @return - a GerepObject
	 * @throws IOException
	 * @throws MethodNotFoundException
	 * @throws ClassNotFoundException
	 */
	public static <T extends GenrepObject> GenrepObject doQuery(String servUrl,METHOD method,FORMAT format,Class<? extends GenrepObject> clazz,KEY key,String query) throws IOException, MethodNotFoundException, ClassNotFoundException{

		String url = prepareUrl(servUrl, method, format, key, query);
		System.out.println("fetching  "+url);
		String result = InternetConnection.sendGETConnection(url);
		GenrepObject genrepObject = JsonMapper.deserializeObject(result, clazz);
		return genrepObject.getInstance();
	}

	/**
	 * Launch the query to Genrep and get
	 * back an GenrepObject
	 * @param servUrl - the Genrep server URL
	 * @param method - the methods supported by Genrep
	 * @param format - the format supported by Genrep
	 * @param clazz - the class you want to get back (map the deserialization of the JSON)
	 * @param key - the key to put in the URL (usually the same than the className)
	 * @param query - the query (e.g : an id,md5,...)
	 * @return - a GerepObject
	 * @throws IOException
	 * @throws MethodNotFoundException
	 * @throws ClassNotFoundException
	 */
	public static  <T extends GenrepObject> GenrepObject doQuery(String servUrl, METHOD method, FORMAT format,
			Class<? extends GenrepObject> clazz, KEY key, int query) throws IOException, MethodNotFoundException, ClassNotFoundException {
		return doQuery(servUrl, method, format, clazz, key, Integer.toString(query));
	}


	public static <T extends GenrepObject> List<? extends GenrepObject> doQueryList(String servUrl, METHOD method, FORMAT format,
			TypeReference<List<GenrepObject>> typeReference, KEY key, String query) throws MethodNotFoundException, IOException {
		String url = prepareUrl(servUrl, method, format, key, query);
		System.out.println("fetching  "+url);
		String result = InternetConnection.sendGETConnection(url);
		List<? extends GenrepObject> genrepObjects = JsonMapper.deserializeArray(result, typeReference);
		List<GenrepObject> newObjects = new ArrayList<GenrepObject>();
		for( GenrepObject o : genrepObjects){
			newObjects.add(o.getInstance());
		}
		return newObjects;
	}




}
