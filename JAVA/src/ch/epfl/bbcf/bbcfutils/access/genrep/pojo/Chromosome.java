package ch.epfl.bbcf.bbcfutils.access.genrep.pojo;

import org.codehaus.jackson.annotate.JsonCreator;


public class Chromosome extends GenrepObject{

	
	@JsonCreator
	public Chromosome() {
	}
	private String name,updated_at,refseq_locus,synonyms,chr_name,created_at;
	private int genome_id,chr_type_id,refseq_version,id,num,length,chromosome_id,assembly_id,gi_number;
	public void setUpdated_at(String updated_at) {
		this.updated_at = updated_at;
	}
	public String getUpdated_at() {
		return updated_at;
	}
	public void setRefseq_locus(String refseq_locus) {
		this.refseq_locus = refseq_locus;
	}
	public String getRefseq_locus() {
		return refseq_locus;
	}
	public void setSynonyms(String synonyms) {
		this.synonyms = synonyms;
	}
	public String getSynonyms() {
		return synonyms;
	}
	public void setChr_name(String chr_name) {
		this.chr_name = chr_name;
	}
	public String getChr_name() {
		return chr_name;
	}
	public void setCreated_at(String created_at) {
		this.created_at = created_at;
	}
	public String getCreated_at() {
		return created_at;
	}
	public void setName(String name) {
		this.name = name;
	}
	public String getName() {
		return name;
	}
	public void setGi_number(int gi_number) {
		this.gi_number = gi_number;
	}
	public int getGi_number() {
		return gi_number;
	}
	public void setChromosome_id(int chromosome_id) {
		this.chromosome_id = chromosome_id;
	}
	public int getChromosome_id() {
		return chromosome_id;
	}
	public void setRefseq_version(int refseq_version) {
		this.refseq_version = refseq_version;
	}
	public int getRefseq_version() {
		return refseq_version;
	}
	public void setChr_type_id(int chr_type_id) {
		this.chr_type_id = chr_type_id;
	}
	public int getChr_type_id() {
		return chr_type_id;
	}
	public void setAssembly_id(int assembly_id) {
		this.assembly_id = assembly_id;
	}
	public int getAssembly_id() {
		return assembly_id;
	}
	public void setGenome_id(int genome_id) {
		this.genome_id = genome_id;
	}
	public int getGenome_id() {
		return genome_id;
	}
	public void setNum(int num) {
		this.num = num;
	}
	public int getNum() {
		return num;
	}
	public void setId(int id) {
		this.id = id;
	}
	public int getId() {
		return id;
	}
	public void setLength(int length) {
		this.length = length;
	}
	public int getLength() {
		return length;
	}
	@Override
	public Chromosome getInstance() {
		return (Chromosome) this.instance;
	}
}
