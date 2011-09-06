package ch.epfl.bbcf.bbcfutils.access.utility;

public class SelectOption {

	private int key;
	private String value;
	public SelectOption(int key_, String value_) {
		this.key = key_;
	}
	public void setKey(int key) {
		this.key = key;
	}
	public int getKey() {
		return key;
	}
	public void setValue(String value) {
		this.value = value;
	}
	public String getValue() {
		return value;
	}
}
