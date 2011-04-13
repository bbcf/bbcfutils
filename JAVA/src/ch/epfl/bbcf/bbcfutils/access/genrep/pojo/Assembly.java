package ch.epfl.bbcf.bbcfutils.access.genrep.pojo;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import ch.epfl.bbcf.bbcfutils.access.genrep.Constants;
import ch.epfl.bbcf.bbcfutils.access.genrep.GenRepAccess;
import ch.epfl.bbcf.bbcfutils.access.genrep.Constants.FORMAT;
import ch.epfl.bbcf.bbcfutils.access.genrep.Constants.KEY;
import ch.epfl.bbcf.bbcfutils.access.genrep.Constants.METHOD;
import ch.epfl.bbcf.bbcfutils.access.genrep.MethodNotFoundException;


public class Assembly extends GenrepObject{

	
	private String name,updated_at,source_name,md5,created_at;
	private boolean bbcf_valid;
	private int nr_assembly_id,genome_id,id,source_id;
	private List<GenrepObject> chromosomes;
	
	public Assembly(){
	}
	public void setCreated_at(String created_at) {
		this.created_at = created_at;
	}
	public String getCreated_at() {
		return created_at;
	}
	public void setUpdated_at(String updated_at) {
		this.updated_at = updated_at;
	}
	public String getUpdated_at() {
		return this.updated_at;
	}
	public void setMd5(String md5) {
		this.md5 = md5;
	}
	public String getMd5() {
		return this.md5;
	}
	public void setName(String name) {
		this.name = name;
	}
	public String getName() {
		return this.name;
	}
	public void setSource_name(String source_name) {
		this.source_name = source_name;
	}
	public String getSource_name() {
		return this.source_name;
	}
	public void setBbcf_valid(boolean bbcf_valid) {
		this.bbcf_valid = bbcf_valid;
	}
	public boolean isBbcf_valid() {
		return this.bbcf_valid;
	}
	public void setId(int id) {
		this.id = id;
	}
	public int getId() {
		return this.id;
	}
	public void setNr_assembly_id(int nr_assembly_id) {
		this.nr_assembly_id = nr_assembly_id;
	}
	public int getNr_assembly_id() {
		return this.nr_assembly_id;
	}
	public void setSource_id(int source_id) {
		this.source_id = source_id;
	}
	public int getSource_id() {
		return this.source_id;
	}
	public void setGenome_id(int genome_id) {
		this.genome_id = genome_id;
	}
	public int getGenome_id() {
		return this.genome_id;
	}
	public void setChromosomes(List<GenrepObject> chromosomes) {
		this.chromosomes = chromosomes;
	}
	public List<Chromosome> getChromosomes() {
		List<Chromosome> chrs = new ArrayList<Chromosome>();
		if(null==chromosomes){
			try {
				Assembly tmp = (Assembly) GenRepAccess.doQuery(
						Constants.URL, METHOD.SHOW, FORMAT.json, Assembly.class,KEY.assemblies, this.id);
				this.chromosomes=tmp.chromosomes;
			} catch (IOException e) {
				e.printStackTrace();
			} catch (MethodNotFoundException e) {
				e.printStackTrace();
			} catch (ClassNotFoundException e) {
				e.printStackTrace();
			}
			
		}
		for(GenrepObject o : chromosomes){
			chrs.add((Chromosome) o.getInstance());
		}
		return chrs;
	}
	@Override
	public Assembly getInstance() {
		return (Assembly) this.instance;
	}
}
