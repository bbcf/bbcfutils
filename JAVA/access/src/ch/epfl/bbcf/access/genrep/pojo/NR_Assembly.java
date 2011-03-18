package ch.epfl.bbcf.access.genrep.pojo;

public class NR_Assembly extends GenrepObject{

	private String name,updated_at,ensembl_release_date,gtf_file_ftp,md5,created_at,cdna_file_ftp;
	private int genome_id,id;
	
	public void setMd5(String md5) {
		this.md5 = md5;
	}
	public String getMd5() {
		return md5;
	}
	public void setName(String name) {
		this.name = name;
	}
	public String getName() {
		return name;
	}
	public void setCreated_at(String created_at) {
		this.created_at = created_at;
	}
	public String getCreated_at() {
		return created_at;
	}
	public void setGtf_file_ftp(String gtf_file_ftp) {
		this.gtf_file_ftp = gtf_file_ftp;
	}
	public String getGtf_file_ftp() {
		return gtf_file_ftp;
	}
	public void setCdna_file_ftp(String cdna_file_ftp) {
		this.cdna_file_ftp = cdna_file_ftp;
	}
	public String getCdna_file_ftp() {
		return cdna_file_ftp;
	}
	public void setUpdated_at(String updated_at) {
		this.updated_at = updated_at;
	}
	public String getUpdated_at() {
		return updated_at;
	}
	public void setEnsembl_release_date(String ensembl_release_date) {
		this.ensembl_release_date = ensembl_release_date;
	}
	public String getEnsembl_release_date() {
		return ensembl_release_date;
	}
	public void setId(int id) {
		this.id = id;
	}
	public int getId() {
		return id;
	}
	public void setGenome_id(int genome_id) {
		this.genome_id = genome_id;
	}
	public int getGenome_id() {
		return genome_id;
	}
	@Override
	public NR_Assembly getInstance() {
		return (NR_Assembly) this.instance;
	}
	
	
}
