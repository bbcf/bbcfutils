package ch.epfl.bbcf.access.genrep;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;

import ch.epfl.bbcf.access.InternetConnection;

public class OrganismsAccess extends GenRepAccess{

	private JSONArray organisms;

	public OrganismsAccess() throws IOException, JSONException{
		String result = InternetConnection.sendGETConnection(GENEREP+ORGANISMS_URL+GET);
		this.organisms = new JSONArray(result);
	}


	/**
	 * get an organism by it's id
	 * @param id - the organism identifier
	 * @return - the JSONObject representing the organism
	 * @throws JSONException
	 */
	public JSONObject getOrganismById(int id) throws JSONException {
		for(int i=0;i<organisms.length();i++){
			JSONObject json_organism = organisms.getJSONObject(i);
			JSONObject organism = json_organism.getJSONObject(ORGANISM_KEY);
			if(json_organism.getInt(ID)==id){
				return organism;
			}
		}
		return no_json_result;
	}




	/**
	 * take the list of organisms in Genrep
	 * @return a list of selectOption (id,name)
	 * @throws JSONException 
	 */
	public List<JSONObject> getOrganisms() throws JSONException {
		List<JSONObject> orgs = new ArrayList<JSONObject>();
		for(int i=0;i<organisms.length();i++){
			JSONObject json_organism = organisms.getJSONObject(i);
			JSONObject organism = json_organism.getJSONObject(ORGANISM_KEY);
			orgs.add(organism);
		}
		return orgs;
	}
}
