package ch.epfl.bbcf.access.genrep;

public interface GenRepKeys {

	public static final String GENEREP = "http://bbcftools.vital-it.ch/genrep/";

	/**
	 * URI
	 */
	public static final String CHROMOSOMES_URL="chromosomes";
	public static final String ASSEMBLIES_URL="assemblies";
	public static final String NR_ASSEMBLIES_URL="nr_assemblies";
	public static final String GENOMES_URL="genomes";
	public static final String ORGANISMS_URL="organisms";
	public static final String REPOSITORY_URL="repo_file";

	public static final String GET=".json";

	/**
	 * Principal keys
	 */
	public final static String GENOME_KEY ="genome";
	public static final String NR_ASSEMBLY_KEY ="nr_assembly";
	public static final String CHROMOSOME_KEY = "chromosomes";
	public final static String ORGANISM_KEY="organism";
	public static final String ASSEMBLY_KEY ="assembly";
	public static final String REPOSITORY_KEY ="repo_file";
	
	/**
	 * General Keys
	 */
	public final static String ID = "id";
	public static final String LENGTH = "length";
	public static final String MD5 = "md5";
	public static final String NAME = "name";
	public static final String SPECIES = "species";
	
	/**
	 * Particulars keys
	 */
	public final static String ORGANISM_ID="organism_id";
	public static final String GENOME_ID ="genome_id";
	public static final String NR_ASSEMBLIES_ID ="nr_assembly_id";
	

	
}
