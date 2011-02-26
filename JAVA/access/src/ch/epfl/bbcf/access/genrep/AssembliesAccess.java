package ch.epfl.bbcf.access.genrep;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;

import ch.epfl.bbcf.access.InternetConnection;

public class AssembliesAccess extends GenRepAccess{

	private JSONArray assemblies;

	public AssembliesAccess() throws IOException, JSONException{
		String result = InternetConnection.sendGETConnection(GENEREP+ASSEMBLIES_URL+GET);
		this.assemblies = new JSONArray(result);
	}



	public JSONObject getGenomeById(int id) throws JSONException {
		for(int i=0;i<assemblies.length();i++){
			JSONObject json_assembly = assemblies.getJSONObject(i);
			JSONObject assembly = json_assembly.getJSONObject(ASSEMBLY_KEY);
			if(json_assembly.getInt(ID)==id){
				return assembly;
			}
		}
		return no_json_result;
	}

	/**
	 * take an array of chromosomes belonging to an assembly
	 * @param id - the nr_assembly identifier
	 * @return
	 * @throws JSONException
	 * @throws IOException
	 */
	public JSONArray getChromosomeListByNR_AssemblyId(int id) throws JSONException, IOException{
		for(int i=0;i<assemblies.length();i++){
			JSONObject json_assembly = assemblies.getJSONObject(i);
			if(json_assembly.getInt(NR_ASSEMBLIES_ID)==id){
				int assemblyId = json_assembly.getInt(ID);
				JSONObject showAssembly = new JSONObject(
						InternetConnection.sendGETConnection(GENEREP+ASSEMBLIES_URL+SHOW(assemblyId))
				);
				JSONObject assembly = showAssembly.getJSONObject(ASSEMBLY_KEY);
				return assembly.getJSONArray(CHROMOSOME_KEY);
			}
		}
		return no_json_array_result;
	}
}