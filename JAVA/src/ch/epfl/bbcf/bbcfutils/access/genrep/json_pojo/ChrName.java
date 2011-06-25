package ch.epfl.bbcf.bbcfutils.access.genrep.json_pojo;

public class ChrName extends GenrepObject{
	
	private String updated_at,value,created_at;
	private int id,chromosome_id,assembly_id;

	

	
	@Override
	public Chromosome getInstance() {
		return (Chromosome) this.instance;
	}




	public void setValue(String value) {
		this.value = value;
	}




	public String getValue() {
		return value;
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
		return updated_at;
	}




	public void setId(int id) {
		this.id = id;
	}




	public int getId() {
		return id;
	}




	public void setAssembly_id(int assembly_id) {
		this.assembly_id = assembly_id;
	}




	public int getAssembly_id() {
		return assembly_id;
	}




	public void setChromosome_id(int chromosome_id) {
		this.chromosome_id = chromosome_id;
	}




	public int getChromosome_id() {
		return chromosome_id;
	}
}
