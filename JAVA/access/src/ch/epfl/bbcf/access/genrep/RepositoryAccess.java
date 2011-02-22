package ch.epfl.bbcf.access.genrep;

import java.io.IOException;

import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;

import ch.epfl.bbcf.access.InternetConnection;

public class RepositoryAccess extends GenRepAccess{

	private JSONArray repositories;

	public RepositoryAccess() throws IOException, JSONException{
		String result = InternetConnection.sendGETConnection(GENEREP+REPOSITORY_URL+GET);
		this.repositories = new JSONArray(result);
	}


	/**
	 * take the URL where to fetch 
	 * the fasta file for the assembly
	 * @param id - the nr_assembly
	 * @return - the url
	 * @throws JSONException
	 * @throws IOException
	 */
	public static String getURLFastaFileByNR_AssemblyId(int id) throws JSONException, IOException{
		NR_AssembliesAccess nr_access = GenRepAccess.getNR_assemblies();
		String md5 = nr_access.getMD5ByNR_AssemblyId(id);
		return GENEREP+"data/nr_assemblies/fasta/"+md5+".tar.gz";
	}

}
