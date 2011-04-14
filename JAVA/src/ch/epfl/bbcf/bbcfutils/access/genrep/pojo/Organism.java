package ch.epfl.bbcf.bbcfutils.access.genrep.pojo;

public class Organism extends GenrepObject{

	private String updated_at,species,created_at;
	private int id,tax_id;
	
	public void setSpecies(String species) {
		this.species = species;
	}
	public String getSpecies() {
		return species;
	}
	public void setUpdated_at(String updated_at) {
		this.updated_at = updated_at;
	}
	public String getUpdated_at() {
		return updated_at;
	}
	public void setCreated_at(String created_at) {
		this.created_at = created_at;
	}
	public String getCreated_at() {
		return created_at;
	}
	public void setTax_id(int tax_id) {
		this.tax_id = tax_id;
	}
	public int getTax_id() {
		return tax_id;
	}
	public void setId(int id) {
		this.id = id;
	}
	public int getId() {
		return id;
	}
	@Override
	public Organism getInstance() {
		return (Organism) this.instance;
	}
	
}
