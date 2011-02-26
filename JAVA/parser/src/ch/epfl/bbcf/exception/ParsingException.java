package ch.epfl.bbcf.exception;

public class ParsingException extends Exception {

	public ParsingException(Exception e) {
		super(e);
	}

	public ParsingException(Exception e, String message, int line) {
		super(e.getLocalizedMessage()+"  "+message+" at line "+line);
	}

	public ParsingException(String message, int line) {
		super(message+" at line "+line);
	}

}
