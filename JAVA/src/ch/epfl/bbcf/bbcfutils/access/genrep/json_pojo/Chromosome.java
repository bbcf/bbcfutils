package ch.epfl.bbcf.bbcfutils.access.genrep.json_pojo;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.codehaus.jackson.annotate.JsonCreator;

import ch.epfl.bbcf.bbcfutils.access.genrep.Constants;
import ch.epfl.bbcf.bbcfutils.access.genrep.GenRepAccess;
import ch.epfl.bbcf.bbcfutils.access.genrep.MethodNotFoundException;
import ch.epfl.bbcf.bbcfutils.access.genrep.Constants.FORMAT;
import ch.epfl.bbcf.bbcfutils.access.genrep.Constants.KEY;
import ch.epfl.bbcf.bbcfutils.access.genrep.Constants.METHOD;


public class Chromosome extends GenrepObject{

	
	@JsonCreator
	public Chromosome() {
	}
	private String name,updated_at,refseq_locus,synonyms,created_at,chr_name;
	private int genome_id,chr_type_id,refseq_version,id,num,length,chromosome_id,assembly_id,gi_number;
	private List<GenrepObject> chr_names;
	
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
	public void setChr_names(List<GenrepObject> chr_names) {
		this.chr_names = chr_names;
	}
	public List<ChrName> getChr_names() {
		List<ChrName> chrNames = new ArrayList<ChrName>();
		for(GenrepObject o : chr_names){
			chrNames.add((ChrName) o.getInstance());
		}
		return chrNames;
	}
	public void setChr_name(String chr_name) {
		this.chr_name = chr_name;
	}
	public String getChr_name() {
		return chr_name;
	}
}
