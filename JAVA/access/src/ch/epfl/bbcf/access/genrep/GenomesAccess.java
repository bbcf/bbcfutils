package ch.epfl.bbcf.access.genrep;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;

import ch.epfl.bbcf.access.InternetConnection;

public class GenomesAccess extends GenRepAccess{

	private JSONArray genomes;

	public GenomesAccess() throws IOException, JSONException{
		String result = InternetConnection.sendGETConnection(GENEREP+GENOMES_URL+GET);
		this.genomes = new JSONArray(result);
	}


	/**
	 * take a genome by it's id
	 * @param id - the genome id
	 * @return - a JSONObject representing the id
	 * @throws JSONException
	 */
	public JSONObject getGenomeById(int id) throws JSONException {
		for(int i=0;i<genomes.length();i++){
			JSONObject json_genome = genomes.getJSONObject(i);
			JSONObject genome = json_genome.getJSONObject(GENOME_KEY);
			if(json_genome.getInt(ID)==id){
				return genome;
			}
		}
		return no_json_result;
	}

	
	/**
	 * take the genomes belonging to an organism
	 * @param id - the organism identifier
	 * @return - a list of JSONObject representing the genomes
	 * @throws JSONException 
	 */
	public List<JSONObject> getGenomesByOrganismId(int id) throws JSONException {
		List<JSONObject> gen = new ArrayList<JSONObject>();
		for(int i=0;i<genomes.length();i++){
			JSONObject json_genome = genomes.getJSONObject(i);
			if(json_genome.getInt(ORGANISM_ID)==id){
				gen.add(json_genome);
			}
		}
		return gen;
	}
	
	
	
	
}
