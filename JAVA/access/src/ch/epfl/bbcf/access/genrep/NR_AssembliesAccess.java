package ch.epfl.bbcf.access.genrep;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;

import ch.epfl.bbcf.access.InternetConnection;


public class NR_AssembliesAccess extends GenRepAccess{

	private JSONArray nr_assemblies;

	/**
	 * initialize the connection to GenRep
	 * and fetch the nr_assemblies.json
	 * @throws JSONException
	 * @throws IOException
	 */
	public NR_AssembliesAccess() throws JSONException, IOException{
		String result = InternetConnection.sendGETConnection(GENEREP+NR_ASSEMBLIES_URL+GET);
		this.nr_assemblies = new JSONArray(result);
	}



	/**
	 * take the JSON nr_assembly by it's id
	 * @param id - the id
	 * @return a JSONObject
	 * @throws JSONException
	 */
	public JSONObject getNr_AssemblyById(int id) throws JSONException{
		for(int i=0;i<nr_assemblies.length();i++){
			JSONObject json_assembly = nr_assemblies.getJSONObject(i);
			JSONObject assembly = json_assembly.getJSONObject(NR_ASSEMBLY_KEY);
			if(assembly.getInt(ID)==id){
				return assembly;
			}
		}
		return no_json_result;
	}

	/**
	 * get all the assemblies belonging to 
	 * an organism
	 * @param id - the organism identifier
	 * @return - a list of JSONObject representing the assemblies
	 * @throws JSONException 
	 * @throws IOException 
	 */
	public List<JSONObject> getNr_assembliesByOrganimId(int id) throws IOException, JSONException{
		GenomesAccess genomeAccess = GenRepAccess.getGenomeAccess();
		List<JSONObject> genomes = genomeAccess.getGenomesByOrganismId(id);
		List<JSONObject> ass = new ArrayList<JSONObject>();
		for(int j=0;j<genomes.size();j++){
			for(int i=0;i<nr_assemblies.length();i++){
				JSONObject gen = genomes.get(j);
				int genId = gen.getInt(ID);
				JSONObject json_assembly = nr_assemblies.getJSONObject(i);
				if(json_assembly.getInt(GENOME_ID)==genId){
					ass.add(json_assembly);
				}
			}
		}
		return ass;

	}

	/**
	 * take the md5 of an nr_assembly with it's id
	 * @param id - the nr_assembly identifier
	 * @return - the md5
	 * @throws JSONException 
	 */
	public String getMD5ByNR_AssemblyId(int id) throws JSONException{
		JSONObject assembly = getNr_AssemblyById(id);
		return assembly.getString(MD5);
	}
	
	
	
	
	
	
	
	
}