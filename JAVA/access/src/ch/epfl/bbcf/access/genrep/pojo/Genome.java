package ch.epfl.bbcf.access.genrep.pojo;

public class Genome extends GenrepObject{

	private String name,updated_at,strain,mutant,created_at;
	private int organism_id,id,tax_id;
	private boolean in_refseq;
	
	public void setName(String name) {
		this.name = name;
	}
	public String getName() {
		return name;
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
	public void setStrain(String strain) {
		this.strain = strain;
	}
	public String getStrain() {
		return strain;
	}
	public void setMutant(String mutant) {
		this.mutant = mutant;
	}
	public String getMutant() {
		return mutant;
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
	public void setOrganism_id(int organism_id) {
		this.organism_id = organism_id;
	}
	public int getOrganism_id() {
		return organism_id;
	}
	public void setIn_refseq(boolean in_refseq) {
		this.in_refseq = in_refseq;
	}
	public boolean isIn_refseq() {
		return in_refseq;
	}
	@Override
	public Genome getInstance() {
		return (Genome) this.instance;
	}
	
}
