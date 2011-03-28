package ch.epfl.bbcf.access.genrep;

import java.util.HashMap;
import java.util.Map;

public class Constants {

	/**
	 * URL of the server
	 */
	public static final String URL ="http://bbcftools.vital-it.ch/genrep";
	
	/**
	 * Different Objects that Genrep support
	 */
	public static final String ASSEMBLY="assembly",
	CHROMOSOME = "chromosome",
	GENOME = "genome",
	ORGANISM = "organism",
	NR_ASSEMBLY = "nr_assembly",
	SOURCE = "source";

	/**
	 * Keys to put in the URL to get
	 * the Genrep Objects
	 */
	public static enum KEY{assemblies,organisms,genomes,nr_assemblies,chromosomes,sources}
	/**
	 * Methods supported by Genrep
	 */
	public static enum METHOD {INDEX,SHOW,ALL};
	/**
	 * Formats supported by Genrep 
	 */
	public static enum FORMAT {json,html,gtf};
	
	/**
	 * List of parameters supported by Genrep
	 * for each methods for the INDEX method
	 */
	public static final Map<KEY,String[]> parameters_supported = buildMapParameters();

	public static final String[] assembliesParameters = {"name","genome_id"},
	repo_files_parameters = {"assembly_id","chromosome_id","file_type_id"},
	chromosomes_parameters = {"genome_id","assembly_id","assembly_name","chromosome_id"};
	

	private static final Map<KEY, String[]> buildMapParameters() {
		Map<KEY,String[]> map = new HashMap<Constants.KEY, String[]>();
		map.put(KEY.assemblies, assembliesParameters);
		map.put(KEY.chromosomes, chromosomes_parameters);
		return null;
	}
}
