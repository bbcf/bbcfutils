package ch.epfl.bbcf.bbcfutils.exception;

public class ConvertToJSONException extends Exception{

	public ConvertToJSONException(String mess) {
		super(mess);
	}

	public ConvertToJSONException(Exception e) {
		super(e);
	}

}
